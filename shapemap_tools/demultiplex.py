import fire
import subprocess as sp
import os
import pandas as pd
import glob
import logging

logger = logging.getLogger()

output_pattern = "{output_path}/{prefix}_{{name}}.fastq"


def basename_without_ext(path):
    return (os.path.basename(path).split(os.extsep))[0]


def run_cutadapt_demultiplex_pairend(
    readfq: str,
    tagfile_path,
    outreadfq=None,
    matefq: str = None,
    outmatefq=None,
    cores=8,
    max_errors=0,
    indels=False,
):
    # Here R1 and  R2 are inverted so the R2 5' primer (which is the 3' primer on the
    # main strand) is used for selecting pairs.
    cmd = [
        "cutadapt",
        "-j",
        str(cores),
        "--action",
        "none",
        "-e",
        str(max_errors),
        "--no-indels" if not indels else "--indels",
        "-g",
        f"file:{tagfile_path}",
        "-o",
        outreadfq,
        "-p",
        outmatefq,
        readfq,
        matefq,
    ]

    sp.run(cmd)


pairs = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def rev_complement(seq):
    rc = ""
    for nuc in seq:
        rc = pairs[nuc] + rc
    return rc


def prepare_tagfile(tagfile, groupcolumn, output_path, max_errors=0, indels=False):
    indels = "indels" if indels else "noindels"
    tagdf = pd.read_csv(tagfile, sep="\t")

    outtag = pd.DataFrame(tagdf[["name"]])

    outtag["tag"] = tagdf["tag"]  # + tagdf["primer"]#.str.slice(0, 4)
    outtag["mate_tag"] = (tagdf["tag"]).apply(rev_complement)  # GCTGGTGTTGCCGTCGA
    outtag["primer"] = tagdf["primer"]
    outtag["revprimer"] = tagdf["primer"].apply(rev_complement)

    if groupcolumn:
        outtag["group"] = tagdf[groupcolumn]

    tagfile_path = os.path.join(output_path, "tags.tsv")
    outtag.to_csv(tagfile_path, sep="\t")  # , index=False, header=False)

    with open(os.path.join(output_path, "tags.fa"), "w") as fd, open(
        os.path.join(output_path, "tags_mate.fa"), "w"
    ) as fdm:
        for rid, row in outtag.iterrows():
            fd.write(
                f">{row['name']}\n{row['tag']}"
                f";e={max_errors};{indels};min_overlap={len(row['tag'])}\n"
            )
            fdm.write(
                f">{row['name']}\n{row['mate_tag']}"
                f";e={max_errors};{indels};min_overlap={len(row['mate_tag'])}\n"
            )

    return tagfile_path[:-4]


def get_fastq_files(folder):
    R1 = glob.glob(f"{folder}/**/*R1_[0-9][0-9][0-9].fastq.gz") + glob.glob(
        f"{folder}/*R1_[0-9][0-9][0-9].fastq.gz"
    )
    R2 = glob.glob(f"{folder}/**/*R2_[0-9][0-9][0-9].fastq.gz") + glob.glob(
        f"{folder}/*R2_[0-9][0-9][0-9].fastq.gz"
    )
    return R1, R2


def get_subpath(folder, fastq):
    subpath = os.path.dirname(fastq.removeprefix(folder.removesuffix("/") + "/"))
    return subpath


def get_output_fastq_path(output_path, subpath, fastq):
    outreadfq = output_pattern.format(
        output_path=os.path.join(output_path, subpath),
        prefix=basename_without_ext(fastq),
    )
    return outreadfq


def makesubdirs(path, subpath):
    if subpath not in ["", "/"]:
        os.makedirs(os.path.join(path, subpath), exist_ok=True)


def demultiplex(
    tagfile,
    folder: str = None,
    R1: [str] = None,
    R2: [str] = None,
    groupcolumn=None,
    output_path="output",
    cores=8,
    max_errors=1,
    indels=False,
):
    os.makedirs(output_path, exist_ok=True)
    tagfile_path = prepare_tagfile(
        tagfile,
        groupcolumn,
        output_path,
        max_errors=max_errors,
        indels=indels,
    )

    if folder:
        if R1 is not None or R2 is not None:
            raise fire.Error("You cannot mix --folder and --R1/--R2 arguments")
        R1, R2 = get_fastq_files(folder)
    if R2:
        assert len(R1) == len(R2)

    if R2:
        logger.info("Using Pairwise demultiplexing")
        for readfq, matefq in zip(R1, R2):

            readfq_subpath = get_subpath(folder, readfq)
            matefq_subpath = get_subpath(folder, matefq)
            makesubdirs(output_path, readfq_subpath)
            makesubdirs(output_path, matefq_subpath)

            outreadfq = get_output_fastq_path(output_path, readfq_subpath, readfq)
            outmatefq = get_output_fastq_path(output_path, matefq_subpath, matefq)

            logger.info("Demultiplexing {os.path.basename(readfq)}")

            run_cutadapt_demultiplex_pairend(
                readfq=readfq,
                matefq=matefq,
                tagfile_path=tagfile_path + "_mate.fa",
                # tagfile_mate_path=tagfile_path + "_mate.fa",
                outreadfq=outreadfq,
                outmatefq=outmatefq,
                cores=cores,
                indels=indels,
                max_errors=max_errors,
            )
            unmapped_readfq = outreadfq.format(name="unknown")
            unmapped_matefq = outmatefq.format(name="unknown")
            outreadfq_rev = outreadfq.format(name="{name}_rev")
            outmatefq_rev = outmatefq.format(name="{name}_rev")
            run_cutadapt_demultiplex_pairend(
                readfq=unmapped_matefq,
                matefq=unmapped_readfq,
                # tagfile_path=tagfile_path + ".fa",
                tagfile_path=tagfile_path + "_mate.fa",
                outreadfq=outmatefq_rev,
                outmatefq=outreadfq_rev,
                cores=cores,
                indels=indels,
                max_errors=max_errors,
            )

        # os.remove(trimmed_readfq)
        # os.remove(trimmed_matefq)


#    else:
#        for readfq in R1:
#            intermed_readfq = os.path.dirname(
#                readfq.removeprefix(folder.removesuffix("/") + "/")
#            )
#            if intermed_readfq != "" and intermed_readfq != "/":
#                os.makedirs(os.path.join(output_path, intermed_readfq), exist_ok=True)
#            trimmed_readfq = os.path.join(
#                output_path, intermed_readfq, "trimmed_" + os.path.basename(readfq)
#            )
#            logger.info("Trimming end {os.path.basename(readfq)}")
#
#            outreadfq = output_pattern.format(
#                output_path=os.path.join(output_path, matefq),
#                prefix=basename_without_ext(readfq),
#            )
#            logger.info("Demultiplexing {os.path.basename(readfq)}")
#            run_cutadapt_demultiplex(
#                readfq=trimmed_readfq,
#                tagfile_path=tagfile_path + ".fa",
#                outreadfq=outreadfq,
#                cores=cores,
#            )
#            os.remove(trimmed_readfq)


def main():
    fire.Fire(demultiplex)


if __name__ == "__main__":
    main()

# adapter3p = (
#    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # NEBNext Adapter
#    + "ATTCCTTTATCGGGGTTTGGGGGGTGGGGGATGATAAAATTGGTGTGGGGGGGG"  # Flowcell sequence ?
# )
#
# nebnext_adapters = [
#    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
#    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
#    # "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
# ]

# def run_cutadapt_demultiplex(
#    readfq: str,
#    tagfile_path,
#    outreadfq=None,
#    matefq: str = None,
#    outmatefq=None,
#    cores=8,
# ):
#    cmd = [
#        "cutadapt",
#        "-j",
#        str(cores),
#        "-e",
#        "1",
#        "--no-indels",
#        "-a",
#        f"file:{tagfile_path}",
#        "-o",
#        outreadfq,
#        "-p",
#        outmatefq,
#        readfq,
#        matefq,
#        "--action",
#        "none",
#    ]
#
#    sp.run(cmd)

# def run_cutadapt_trimming(
#    readfq: str, matefq: str = None, outreadfq=None, outmatefq=None, cores=8
# ):
#    cmd = [
#        "cutadapt",
#        "-j",
#        str(cores),
#        "-a",
#        nebnext_adapters[0],
#        "-A",
#        nebnext_adapters[1],
#        "-o",
#        outreadfq,
#        "-p",
#        outmatefq,
#        readfq,
#        matefq,
#    ]
#
#    sp.run(cmd)

# def run_bbduk(readfq: str, matefq: str = None, outreadfq=None, outmatefq=None):
#    cmd = [
#        "bbduk.sh",
#        f"in1={readfq}",
#        f"out1={outreadfq}",
#        "k=21",
#        "mink=8",
#        "ktrim=r",
#        "ref=truseq",
#        # "ref=adapters/TruSeq3-PE-2.fa",
#        "hdist=2",
#        "tpe",
#        "tbo",
#    ]
#
#    if matefq is not None:
#        cmd.extend([f"in2={matefq}", f"out2={outmatefq}"])
#
#    sp.run(cmd)
