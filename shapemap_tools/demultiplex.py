import fire
import subprocess as sp
import os
import pandas as pd
import glob
import logging

logger = logging.getLogger()

output_pattern = "{output_path}/{prefix}_{{name}}.fastq"

adapter3p = (
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # NEBNext Adapter
    + "ATTCCTTTATCGGGGTTTGGGGGGTGGGGGATGATAAAATTGGTGTGGGGGGGG"  # Flowcell sequence ?
)

nebnext_adapters = [
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    # "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
]


def basename_without_ext(path):
    return (os.path.basename(path).split(os.extsep))[0]


def run_cutadapt_trimming(
    readfq: str, matefq: str = None, outreadfq=None, outmatefq=None, cores=8
):
    cmd = [
        "cutadapt",
        "-j",
        str(cores),
        "-a",
        nebnext_adapters[0],
        "-A",
        nebnext_adapters[1],
        "-o",
        outreadfq,
        "-p",
        outmatefq,
        readfq,
        matefq,
    ]

    sp.run(cmd)

print("LOLO")

def run_cutadapt_demultiplex_with_mate(
    readfq: str,
    tagfile_path,
    tagfile_mate_path,
    outreadfq=None,
    matefq: str = None,
    outmatefq=None,
    cores=8,
):
    cmd = [
        "cutadapt",
        "-j",
        str(cores),
        "-e",
        "2",
        "--no-indels",
        "-a",
        f"file:{tagfile_path}",
        "-G",
        f"file:{tagfile_mate_path}",
        "-o",
        outreadfq,
        "-p",
        outmatefq,
        readfq,
        matefq,
    ]

    sp.run(cmd)


def run_cutadapt_demultiplex(
    readfq: str,
    tagfile_path,
    outreadfq=None,
    matefq: str = None,
    outmatefq=None,
    cores=8,
):
    cmd = [
        "cutadapt",
        "-j",
        str(cores),
        "-e",
        "2",
        "--no-indels",
        "-a",
        f"file:{tagfile_path}",
        "-o",
        outreadfq,
        "-p",
        outmatefq,
        readfq,
        matefq,
	"--action",
	"none"
    ]

    sp.run(cmd)


def run_bbduk(readfq: str, matefq: str = None, outreadfq=None, outmatefq=None):
    cmd = [
        "bbduk.sh",
        f"in1={readfq}",
        f"out1={outreadfq}",
        "k=21",
        "mink=8",
        "ktrim=r",
        "ref=truseq",
        # "ref=adapters/TruSeq3-PE-2.fa",
        "hdist=2",
        "tpe",
        "tbo",
    ]

    if matefq is not None:
        cmd.extend([f"in2={matefq}", f"out2={outmatefq}"])

    sp.run(cmd)


pairs = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def rev_complement(seq):
    rc = ""
    for nuc in seq:
        rc = pairs[nuc] + rc
    return rc


def prepare_tagfile(tagfile, groupcolumn, output_path):
    tagdf = pd.read_csv(tagfile, sep="\t")

    outtag = pd.DataFrame(tagdf[["name"]])

    outtag["tag"] = tagdf["tag"]# + tagdf["primer"]#.str.slice(0, 4)
    outtag["mate_tag"] = (tagdf["tag"] + tagdf["primer"]).apply(rev_complement, axis=1) #.str.slice(0, 7)

    if groupcolumn:
        outtag["group"] = tagdf[groupcolumn]

    tagfile_path = os.path.join(output_path, "tags.tsv")
    outtag.to_csv(tagfile_path, sep="\t")#, index=False, header=False)

    with open(os.path.join(output_path, "tags.fa"), "w") as fd, open(os.path.join(output_path, "tags_mate.fa"), "w") as fdm:
        for rid, row in outtag.iterrows():
            fd.write(f">{row['name']}\n{row['tag']}\n")
            fdm.write(f">{row['name']}\n{row['mate_tag']}\n")

    return tagfile_path[:-4]


def demultiplex(
    tagfile,
    folder: str = None,
    R1: [str] = None,
    R2: [str] = None,
    groupcolumn=None,
    output_path="output",
    cores=8,
):
    os.makedirs(output_path, exist_ok=True)
    tagfile_path = prepare_tagfile(tagfile, groupcolumn, output_path)

    if folder:
        if R1 is not None or R2 is not None:
            raise fire.Error("You cannot mix --folder and --R1/--R2 arguments")
        R1 = glob.glob(f"{folder}/**/*R1_[0-9][0-9][0-9].fastq.gz") + glob.glob(
            f"{folder}/*R1_[0-9][0-9][0-9].fastq.gz"
        )
        R2 = glob.glob(f"{folder}/**/*R2_[0-9][0-9][0-9].fastq.gz") + glob.glob(
            f"{folder}/*R2_[0-9][0-9][0-9].fastq.gz"
        )

    if R2:
        assert len(R1) == len(R2)

    if R2:
        logger.info("Using Pairwise demultiplexing")
        for readfq, matefq in zip(R1, R2):
            intermed_readfq = os.path.dirname(
                readfq.removeprefix(folder.removesuffix("/") + "/")
            )
            intermed_matefq = os.path.dirname(
                matefq.removeprefix(folder.removesuffix("/") + "/")
            )
            if intermed_readfq != "" and intermed_readfq != "/":
                os.makedirs(os.path.join(output_path, intermed_readfq), exist_ok=True)
            if intermed_matefq != "" and intermed_matefq != "/":
                os.makedirs(os.path.join(output_path, intermed_matefq), exist_ok=True)

	#            trimmed_readfq = os.path.join(
	#                output_path, intermed_readfq, "trimmed_" + os.path.basename(readfq)
	#            )
	#
	#            trimmed_matefq = os.path.join(
	#                output_path, intermed_matefq, "trimmed_" + os.path.basename(matefq)
	#            )
            outreadfq = output_pattern.format(
                output_path=os.path.join(output_path, intermed_readfq),
                prefix=basename_without_ext(readfq),
            )
            outmatefq = output_pattern.format(
                output_path=os.path.join(output_path, intermed_matefq),
                prefix=basename_without_ext(matefq),
            )

           # logger.info("Trimming end {os.path.basename(readfq)}")
           # run_cutadapt_trimming(
           #     readfq,
           #     matefq,
           #     outreadfq=trimmed_readfq,
           #     outmatefq=trimmed_matefq,
           #     cores=cores,
           # )
            logger.info("Demultiplexing {os.path.basename(readfq)}")

            run_cutadapt_demultiplex_with_mate(
                readfq=readfq,
                matefq=matefq,
                tagfile_path=tagfile_path + ".fa",
                tagfile_mate_path=tagfile_path + "_mate.fa",
                outreadfq=outreadfq,
                outmatefq=outmatefq,
                cores=cores,
            )

           # os.remove(trimmed_readfq)
           # os.remove(trimmed_matefq)
    else:
        for readfq in R1:
            intermed_readfq = os.path.dirname(
                readfq.removeprefix(folder.removesuffix("/") + "/")
            )
            if intermed_readfq != "" and intermed_readfq != "/":
                os.makedirs(os.path.join(output_path, intermed_readfq), exist_ok=True)
            trimmed_readfq = os.path.join(
                output_path, intermed_readfq, "trimmed_" + os.path.basename(readfq)
            )
            logger.info("Trimming end {os.path.basename(readfq)}")
            run_cutadapt_trimming(readfq, outreadfq=trimmed_readfq, cores=cores)

            outreadfq = output_pattern.format(
                output_path=os.path.join(output_path, matefq),
                prefix=basename_without_ext(readfq),
            )
            logger.info("Demultiplexing {os.path.basename(readfq)}")
            run_cutadapt_demultiplex(
                readfq=trimmed_readfq,
                tagfile_path=tagfile_path + ".fa",
                outreadfq=outreadfq,
                cores=cores,
            )
            os.remove(trimmed_readfq)


def main():
    fire.Fire(demultiplex)


if __name__ == "__main__":
    main()
