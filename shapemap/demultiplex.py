import fire
import subprocess as sp
import os
import pandas as pd
import glob

output_pattern = "{output_path}/{prefix}_%.fastq.gz"


# TODO correct path in bbduk
def basename_without_ext(path):
    return os.path.splitext(os.path.basename(path))[0]


def run_bbduk(readfq: str, matefq: str = None, outreadfq=None, outmatefq=None):
    cmd = [
        "bbduk.sh",
        f"in1={readfq}",
        f"out1={outreadfq}",
        "k=25",
        "mink=8",
        "ktrim=r",
        "ref=truseq.fa.gz",
        "hdist=1",
    ]

    if matefq is not None:
        cmd.extend([f"in2={matefq}", f"out2={outmatefq}"])

    sp.run(cmd)


def run_fastq_multx(readfq, outreadfq, tagfile_path, matefq=None, outmatefq=None):

    fastqs = [readfq]
    if matefq is not None:
        fastqs.extend([matefq, "-o", outreadfq, "-o", outmatefq])
    else:
        fastqs.extend(["-o", outreadfq])

    cmd = [
        "fastq-multx",
        # "-e",
        "-b",
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
    tagdf = pd.read_csv(tagfile, sep="\t")

    outtag = pd.DataFrame(tagdf[["name"]])

    outtag["tag"] = tagdf["tag"] + tagdf["primer"]

    if groupcolumn:
        outtag["group"] = tagdf[groupcolumn]

    tagfile_path = os.path.join(output_path, "tags.tsv")
    outtag.to_csv(tagfile_path, sep="\t", index=False, header=False)
    return tagfile_path


def demultiplex(
    tagfile,
    folder: str = None,
    R1: [str] = None,
    R2: [str] = None,
    groupcolumn=None,
    output_path="output",
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
        for readfq, matefq in zip(R1, R2):
            trimmed_readfq = os.path.join(
                os.path.dirname(readfq), "trimmed_" + os.path.basename(readfq)
            )

            trimmed_matefq = os.path.join(
                os.path.dirname(matefq), "trimmed_" + os.path.basename(matefq)
            )
            run_bbduk(
                readfq, matefq, outreadfq=trimmed_readfq, outmatefq=trimmed_matefq
            )

            outreadfq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(readfq)
            )
            outmatefq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(matefq)
            )
            run_fastq_multx(
                readfq=trimmed_readfq,
                matefq=trimmed_matefq,
                tagfile_path=tagfile_path,
                outreadfq=outreadfq,
                outmatefq=outmatefq,
            )
    else:
        for readfq in R1:
            trimmed_readfq = os.path.join(
                os.path.dirname(readfq), "trimmed_" + os.path.basename(readfq)
            )
            run_bbduk(readfq, outreadfq=trimmed_readfq)
            outreadfq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(readfq)
            )
            run_fastq_multx(
                readfq=trimmed_readfq, tagfile_path=tagfile_path, outreadfq=outreadfq
            )


def main():
    fire.Fire(demultiplex)


if __name__ == "__main__":
    main()
