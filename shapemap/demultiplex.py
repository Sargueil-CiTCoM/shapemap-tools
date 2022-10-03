import fire
import subprocess as sp
import os
import pandas as pd
import glob
import logging

logger = logging.getLogger()

output_pattern = "{output_path}/{prefix}_{{name}}.fastq.gz"


# TODO correct path in bbduk
def basename_without_ext(path):
    return os.path.splitext(os.path.basename(path))[0]


def run_cutadapt_trimming(
    readfq: str, matefq: str = None, outreadfq=None, outmatefq=None, cores=8
):
    cmd = [
        "cutadapt",
        "-j",
        str(cores),
        "-a",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "-A",
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        # "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
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
        #        "-e",
        #        "0.15",
        #        "--no-indels",
        "-a",
        f"file:{tagfile_path}",
        "-o",
        outreadfq,
        "-p",
        outmatefq,
        readfq,
        matefq,
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


def run_trimmomatic(
    readfq: str,
    matefq: str = None,
    outreadfq=None,
    outreadfq_unpaired=None,
    outmatefq=None,
    outmatefq_unpaired=None,
):
    if matefq:
        cmd = [
            "trimmomatic",
            "PE",
            "-phred33",
            readfq,
            matefq,
            outreadfq,
            outreadfq_unpaired,
            outmatefq,
            outmatefq_unpaired,
            "ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:4:20:10",
            #            "LEADING:3",
            #            "TRAILING:3",
            # "SLIDINGWINDOW:4:15",
            "MINLEN:36",
        ]
    else:
        cmd = [
            "trimmomatic",
            "SE",
            "-phred33",
            readfq,
            outreadfq,
            outreadfq_unpaired,
            "ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36",
        ]

    sp.run(cmd)


#


def run_fastq_multx(readfq, outreadfq, tagfile_path, matefq=None, outmatefq=None):

    fastqs = [readfq]
    if matefq is not None:
        fastqs.extend([matefq, "-o", outreadfq, "-o", outmatefq])
    else:
        fastqs.extend(["-o", outreadfq])

    cmd = [
        "fastq-multx",
        "-x",
        "-e",
        # "-b",
        "-m1",
        "-B",
        tagfile_path,
    ] + fastqs

    sp.run(cmd)
    # try:
    #    sp.run(cmd)
    # except Exception:
    #    pass


def prepare_tagfile(tagfile, groupcolumn, output_path):
    adapter3p = (
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # NEBNext Adapter
        #+ "ATTCCTTTATCGGGGTTTGGGGGGTGGGGGATGATAAAATTGGTGTGGGGGGGG" # Flowcell sequence ?
    )
    tagdf = pd.read_csv(tagfile, sep="\t")

    outtag = pd.DataFrame(tagdf[["name"]])

    outtag["tag"] = tagdf["tag"] + tagdf["primer"] + adapter3p + "x"

    if groupcolumn:
        outtag["group"] = tagdf[groupcolumn]

    tagfile_path = os.path.join(output_path, "tags.tsv")
    outtag.to_csv(tagfile_path, sep="\t", index=False, header=False)

    with open(os.path.join(output_path, "tags.fa"), "w") as fd:
        for rid, row in outtag.iterrows():
            fd.write(f">{row['name']}\n{row['tag']}\n")

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
        R1 = glob.glob(f"{folder}/*R1_[0-9][0-9][0-9].fastq.gz")
        R2 = glob.glob(f"{folder}/*R2_[0-9][0-9][0-9].fastq.gz")

    if R2:
        assert len(R1) == len(R2)

    if R2:
        logger.info("Using Pairwise demultiplexing")
        for readfq, matefq in zip(R1, R2):
            trimmed_readfq = os.path.join(
                output_path, "trimmed_" + os.path.basename(readfq)
            )

            trimmed_matefq = os.path.join(
                output_path, "trimmed_" + os.path.basename(matefq)
            )
            # trimmed_readfq_unpaired = os.path.join(
            #    output_path, "trimmed_unpaired_" + os.path.basename(readfq)
            # )

            # trimmed_matefq_unpaired = os.path.join(
            #    output_path, "trimmed_unpaired_" + os.path.basename(matefq)
            # )
            logger.info("Trimming end {os.path.basename(readfq)}")
            run_cutadapt_trimming(
                readfq,
                matefq,
                outreadfq=trimmed_readfq,
                outmatefq=trimmed_matefq,
                cores=cores,
            )
            # run_bbduk(
            #    readfq, matefq, outreadfq=trimmed_readfq, outmatefq=trimmed_matefq
            # )
            # run_trimmomatic(
            #    readfq,
            #    matefq,
            #    outreadfq=trimmed_readfq,
            #    outreadfq_unpaired=trimmed_readfq_unpaired,
            #    outmatefq=trimmed_matefq,
            #    outmatefq_unpaired=trimmed_matefq_unpaired,
            # )

            outreadfq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(readfq)
            )
            outmatefq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(matefq)
            )
            logger.info("Demultiplexing {os.path.basename(readfq)}")

            run_cutadapt_demultiplex(
                readfq=trimmed_readfq,
                matefq=trimmed_matefq,
                tagfile_path=tagfile_path + ".fa",
                outreadfq=outreadfq,
                outmatefq=outmatefq,
                cores=cores,
            )
            # run_fastq_multx(
            #    readfq=trimmed_readfq,
            #    matefq=trimmed_matefq,
            #    tagfile_path=tagfile_path,
            #    outreadfq=outreadfq,
            #    outmatefq=outmatefq,
            # )
    else:
        for readfq in R1:
            trimmed_readfq = os.path.join(
                output_path, "trimmed_" + os.path.basename(readfq)
            )
            # trimmed_readfq_unpaired = os.path.join(
            #    output_path, "trimmed_unpaired" + os.path.basename(readfq)
            # )
            logger.info("Trimming end {os.path.basename(readfq)}")
            # run_trimmomatic(
            #    readfq,
            #    outreadfq=trimmed_readfq,
            #    outreadfq_unpaired=trimmed_readfq_unpaired,
            # )
            run_cutadapt_trimming(readfq, outreadfq=trimmed_readfq, cores=cores)

            # run_bbduk(readfq, outreadfq=trimmed_readfq)
            outreadfq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(readfq)
            )
            logger.info("Demultiplexing {os.path.basename(readfq)}")
            run_cutadapt_demultiplex(
                readfq=trimmed_readfq,
                tagfile_path=tagfile_path + ".fa",
                outreadfq=outreadfq,
                cores=cores,
            )

            # run_fastq_multx(
            #    readfq=trimmed_readfq, tagfile_path=tagfile_path, outreadfq=outreadfq
            # )


def main():
    fire.Fire(demultiplex)


if __name__ == "__main__":
    main()
