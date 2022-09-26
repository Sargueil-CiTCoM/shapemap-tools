import fire
import subprocess as sp
import os
import pandas as pd
import glob

output_pattern = "{output_path}/{prefix}_%.fastq.gz"


def basename_without_ext(path):

    return os.path.splitext(os.path.basename(path))[0]


def run_fastq_multx(readfq, outreadfq, tagfile_path, matefq=None, outmatefq=None):

    fastqs = [readfq]
    if matefq is not None:
        fastqs.extend([matefq, "-o", outreadfq, "-o", outmatefq])
    else:
        fastqs.extend(["-o", outreadfq])

    cmd = [
        "fastq_multx",
        "-e",
        "-m1",
        "-B",
        tagfile_path,
    ] + fastqs

    try:
        sp.run(cmd)
    except Exception:
        pass


def prepare_tagfile(tagfile, groupcolumn, output_path):
    tagdf = pd.read_csv(tagfile, sep="\t")

    outtag = tagdf[["name"]]
    outtag["tag"] = tagdf["tag"] + tagdf["primer"]

    if groupcolumn:
        outtag["group"] = tagdf[groupcolumn]

    tagfile_path = os.path.join(output_path, "tags.tsv")
    outtag.write_csv(tagfile_path, sep="\t")
    return tagfile_path


def demultiplex(
    folder: str = None,
    R1: [str] = None,
    R2: [str] = None,
    tagfile="tags.tsv",
    groupcolumn=None,
    output_path="output",
):
    tagfile_path = prepare_tagfile(tagfile, groupcolumn, output_path)

    if folder:
        if R1 is not None or R2 is not None:
            raise fire.Error("You cannot mix --folder and --R1/--R2 arguments")
        R1 = glob.glob(f"{folder}/*R1_[0-9][0-9][0-9].fastq.gz")
        R2 = glob.glob(f"{folder}/*R2_[0-9][0-9][0-9].fastq.gz")
    if R2:
        assert len(R1) == len(R2)
        for readfq, matefq in zip(R1, R2):
            outreadfq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(readfq)
            )
            outmatefq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(matefq)
            )

            run_fastq_multx(
                readfq=readfq,
                matefq=matefq,
                tagfile_path=tagfile_path,
                outreadfq=outreadfq,
                outmatefq=outmatefq,
            )
    else:
        for readfq in R1:
            outreadfq = output_pattern.format(
                output_path=output_path, prefix=basename_without_ext(readfq)
            )
            run_fastq_multx(
                readfq=readfq, tagfile_path=tagfile_path, outreadfq=outreadfq
            )


def main():
    fire.Fire(demultiplex)


if __name__ == "__main__":
    main()
