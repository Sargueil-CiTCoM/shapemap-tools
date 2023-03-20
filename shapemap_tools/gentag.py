#!/usr/bin/env python3

from dataclasses import dataclass
import pandas as pd
import numpy as np
import skbio as sb
import subprocess
import logging
import copy
import seaborn as sns
import matplotlib.pyplot as plt
import fire

sns.set_theme()
pd.set_option("display.max_rows", None)  # or 1000

RNAFOLD_TMP = "tmp/"

logging.basicConfig(level=logging.INFO)


def gen_sequence(size: int = 7):
    rng = np.random.default_rng(42)
    assert size > 0
    base = ["A", "U", "G", "C"]

    seq = "".join([base[rng.integers(0, 3)] for i in range(size)])

    return seq


def load_dna_seq(file):
    sequences = []
    for seq in sb.io.read(file, format="fasta"):
        sequences.append(sb.DNA(seq))
    return sequences


def dist_check(sequences, mindist=5):
    """
    Verify sequence similarity using distance
    Does not take into account shifts and gaps. use with caution.
    """
    count = 0
    ident = 0
    for seq1 in sequences:
        for seq2 in sequences:
            if seq1 != seq2:
                msize = min(len(seq1), len(seq2))
                dist = seq1[:msize].distance(seq2[:msize])
                if dist * msize < mindist:  # and abs(len(seq2) - len(seq1)) < 5:
                    count += 1
                    if dist == 0.0:
                        ident += 1
                        logging.info("IDENTIQUE")
                    logging.info(f"score: {dist} - dist:{dist * msize}")
                    logging.info(
                        f'{len(seq1)}\n>{seq1.metadata["id"]}\n{seq2}\n'
                        f'{len(seq2)}\n>{seq1.metadata["id"]}\n'
                        f"{seq2}\n need a tag\n"
                    )
    return (count, ident)


def eliminate_duplicates(sequences):
    uniq_seq = []
    for iseq1, seq1 in enumerate(sequences):
        uniq = True
        for seq2 in sequences[iseq1 + 1 :]:
            msize = min(len(seq1), len(seq2))
            dist = seq1[:msize].distance(seq2[:msize])
            if dist == 0.0:
                uniq = False
        if uniq:
            uniq_seq.append(seq1)
    return uniq_seq


def rnafold_partition_matrice(
    seq,
    shape_file=None,
    identifier="seq",
    path="./",
):
    shapeflags = ""
    if shape_file:
        shapeflags = f"--shape={shape_file}"

    identifier = identifier.replace(":", "_")
    identifier = identifier.replace("/", "_")
    cmd = (
        f'cd {path} && echo "> {identifier}\n{seq}" | '
        f"RNAfold -p --MEA {shapeflags} > {identifier}.rnafold.out"
    )
    logging.info(cmd)
    process = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    if process.stdout is not None and process.stdout.decode("utf-8") != "":
        logging.info(process.stdout.decode("utf-8"))
    if process.stderr is not None and process.stderr.decode("utf-8") != "":
        logging.error(process.stderr.decode("utf-8"))
    with open(f"{path}{identifier}_dp.ps") as file:
        dotplot = []
        line = file.readline()
        while line and line != "%start of base pair probability data\n":
            line = file.readline()
        line = file.readline()
        while line and line != "showpage\n":
            splt = line[:-1].split(" ")
            ts = [float(ln) for ln in splt[:-1]]
            ts.append(splt[-1])
            dotplot.append(ts)
            line = file.readline()
        coo = pd.DataFrame(dotplot, columns=["x", "y", "proba", "box"])
        return coo
    return None


# def fold_seq(curseq):
#    cmd = f"echo {curseq} | RNAfold -p --MEA"
#    process =subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#    with open('dot.ps') as file:
#        dotplot = []
#
#        line = file.readline()
#        while line and line != "%start of base pair probability data\n":
#            line = file.readline()
#        line = file.readline()
#        while line and line != "showpage\n":
#            splt = line[:-1].split(" ")
#            ts = [float(l) for l in splt[:-1]]
#            ts.append(splt[-1])
#            dotplot.append(ts)
#            line = file.readline()
#
#        coo = pd.DataFrame(dotplot, columns=['x', 'y', 'proba', 'box'])
#        os.delete('dot.ps')
#        return coo
#    os.delete('dot.ps')
#    return None


def dotplot_mat(pds):
    mx = max(int(np.max(pds["x"])), int(np.max(pds["y"])))
    mat = np.zeros((mx, mx))

    for idx, row in pds.iterrows():
        if row["box"] == "ubox":
            mat[int(row["x"] - 1), int(row["y"] - 1)] = row["proba"]
        else:
            mat[int(row["y"] - 1), int(row["x"] - 1)] = row["proba"]
    return mat


def plot_dotplot(dotplot, taggeddotplot):
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 12))
    sns.heatmap(dotplot, ax=ax1, cmap="YlGnBu")
    sns.heatmap(taggeddotplot, ax=ax2, cmap="YlGnBu")
    plt.show()


def accept_tagging(
    orig_seq: sb.RNA,
    tagged_seq: sb.RNA,
    taggeddotplot: pd.DataFrame,
    threshold: float = 0.30,
):
    dpsize = len(orig_seq)

    over_thresdf = taggeddotplot[
        (taggeddotplot["y"] > (dpsize + 7))
        & (taggeddotplot["x"] < (dpsize - 7))
        & (taggeddotplot["proba"] > threshold)
    ]
    if over_thresdf.shape[0] > 0:
        return False, over_thresdf
    else:
        return True, None


@dataclass
class TaggedSequence:
    """Class for keeping track of an item in inventory."""

    sequence: sb.RNA
    tagged_sequence: sb.RNA
    primer: sb.RNA
    tag: sb.RNA
    name: str


def check_tagging(seq: sb.RNA, tagged_seq: sb.RNA, threshold: float = 0.30):
    # if seq:
    #    spares_dotplot = rnafold_partition_matrice(seq, path=RNAFOLD_TMP, identifier=seq.metadata['id'])
    #    #dotplot = dotplot_mat(spares_dotplot)
    #
    spares_dotplotprimer = rnafold_partition_matrice(
        tagged_seq, identifier=seq.metadata["id"], path=RNAFOLD_TMP
    )
    # tagged_dotplot = dotplot_mat(spares_dotplotprimer)

    accepted, maxproba = accept_tagging(
        seq, tagged_seq, spares_dotplotprimer, threshold
    )

    return accepted
    # plot_dotplot(dotplot, tagged_dotplot)
    # print(accepted)
    # print(maxproba)


def tag_sequences(seqs, tags, primers, name_prefix="APSAMN", threshold=0.30):
    tagged_seqs = []
    rejected_seqs = []
    avail_tags = copy.deepcopy(tags)
    taken_tags = []

    for sid, seq in enumerate(seqs):
        accepted_seq = None
        logging.info(f"Starting seq n {sid}")
        for pid, primer in enumerate(primers):
            for tagid, tag in enumerate(avail_tags):
                logging.info(f"Testing tag '{tag}' + primer '{primer}'")
                taggedseq = sb.RNA.concat([seq, tag, primer])
                # spares_dotplot = rnafold_partition_matrice(taggedseq, )
                # dotplot = dotplot_mat(spares_dotplot)
                if check_tagging(seq, taggedseq, threshold):
                    accepted_seq = TaggedSequence(
                        seq, taggedseq, primer, tag, name_prefix + str(sid).zfill(3)
                    )
                    del avail_tags[tagid]
                    taken_tags.append(tag)
                    break
            if accepted_seq:
                break
        if accepted_seq:
            logging.info(f"Accepted sequence {sid}")
            tagged_seqs.append(accepted_seq)
        else:
            logging.warning(f"Rejected sequence {sid}")
            rejected_seqs.append(seq)
        accepted_seq = None

    return tagged_seqs, rejected_seqs, avail_tags, taken_tags


def gentag(
    sequences="seq_ribo_concat.fasta",
    primers="primers.fa",
    tags="tags.list",
    prefix="APSAMN",
    output="output.tsv",
    threshold=0.30,
    remainingtags=None,
):

    rnaseqs = [seq.transcribe() for seq in load_dna_seq(sequences)]
    primers = [seq.transcribe() for seq in load_dna_seq(primers)]
    tagseqs = [seq.transcribe() for seq in load_dna_seq(tags)]

    rnaseqs = eliminate_duplicates(rnaseqs)

    tagged_seqs, rejected_seqs, avail_tags, taken_tags = tag_sequences(
        rnaseqs, tagseqs, primers, prefix, threshold=threshold
    )

    all_tagged_seqs = copy.deepcopy(tagged_seqs)

    tdescs = []
    tids = []
    tnames = []
    tsequences = []
    ttags = []
    tprimers = []
    tsequences_tags_primers = []
    tt7_sequences_tags_primers = []
    t7 = sb.DNA("CGGCGAATCTAATACGACTCACTATAGG").transcribe()

    for idx, seq in enumerate(all_tagged_seqs):
        tdescs.append(seq.sequence.metadata["description"])
        tids.append(seq.sequence.metadata["id"])
        tnames.append(seq.name)
        tsequences.append(str(seq.sequence.reverse_transcribe()))
        ttags.append(seq.tag.reverse_transcribe())
        tprimers.append(seq.primer.reverse_transcribe())
        tsequences_tags_primers.append(str(seq.tagged_sequence.reverse_transcribe()))
        tt7_sequences_tags_primers.append(
            sb.RNA.concat([t7, seq.tagged_sequence]).reverse_transcribe()
        )

    sequence_table = pd.DataFrame(
        {
            "description": tdescs,
            "id": tids,
            "name": tnames,
            "sequence": tsequences,
            "tag": ttags,
            "primer": tprimers,
            "sequence_tag_primer": tsequences_tags_primers,
            "t7_sequence_tag_primer": tt7_sequences_tags_primers,
        }
    )

    sequence_table.to_csv(output, sep="\t")

    if remainingtags:
        with open(remainingtags, "w") as rtags:
            for seq in avail_tags:
                seq.write(rtags)


def main_wrapper():
    fire.Fire(gentag)


if __name__ == "__main__":
    main_wrapper()
