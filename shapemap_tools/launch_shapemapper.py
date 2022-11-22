import fire
import subprocess as sp
import pandas as pd
import ruamel.yaml as yaml
import glob
import copy
from . import fasta
import os
from tqdm import tqdm

shapemapper_path = "/data/fxlyonnet/shapemapper-2.1.5/shapemapper"

YAML = yaml.YAML()

fastq_input_pattern = "{input_path}/{cond}/*_{read}_*.fastq*"
fastq_input_splitted_pattern = "{input_path}/{cond}/*_{read}_*_{sequence}.fastq*"


def run_shapemapper(
    reference: str,
    output_dir: str,
    name: str,
    input_fastqs,
    cores: int = 8,
    indiv_norm: bool = True,
    min_depth: int = 5000,
    log: str = "log.log",
    overwrite=True,
    extra_args: [str] = [],
):
    cmd = [
        shapemapper_path,
        "--name",
        name,
        "--nproc",
        str(cores),
        "--min-depth",
        str(min_depth),
        # "--amplicon",
        "--target",
        reference,
        "--out",
        output_dir,
        "--log",
        log,
    ]
    if indiv_norm:
        cmd.append("--indiv-norm")
    if overwrite:
        cmd.append("--overwrite")

    for condtype in input_fastqs.keys():
        cmd.append(f"--{condtype}")
        for readtype, fastqs in input_fastqs[condtype].items():
            cmd.append(f"--{readtype}")
            cmd.extend(fastqs)

    # print(" ".join(cmd))
    try:
        sp.run(cmd, capture_output=True, text=True)
    except:
        print(f"Error launching shapemapper for {name}")


def fasta_from_tsv(
    tsv_path,
    fasta_path,
    name_col="name",
    seq_col="sequence",
    suffix_cols=["tag", "primer"],
    prefix_cols=[],
):
    os.makedirs(os.path.dirname(fasta_path), exist_ok=True)
    tsvdf = pd.read_csv(tsv_path, sep="\t")
    with open(fasta_path, "w") as fd:
        for rid, row in tsvdf.iterrows():
            name = row[name_col].strip()
            fd.write(f">{name}\n")
            for col in prefix_cols:
                fd.write(row[col].lower())
            fd.write(row[seq_col])
            for col in suffix_cols:
                fd.write(row[col].lower())
            fd.write("\n")


def gen_splitted_ref(ref_path, output_path) -> dict:
    os.makedirs(output_path, exist_ok=True)

    splitted = {}
    ref_path_prefix = os.path.basename(ref_path).split(os.extsep)[0]
    for name, seq in fasta.fasta_iter(ref_path):
        cur_out_path = os.path.join(output_path, f"{ref_path_prefix}_{name}.fasta")
        with open(cur_out_path, "w") as fd:
            fd.write(f">{name}\n")
            fd.write(seq)
        splitted[name] = cur_out_path

    return splitted


def prepare_launch(config, samples):
    runs = {}
    splitted_refs = {}
    afastq = {
        "modified": {"R1": [], "R2": []},
        "untreated": {"R1": [], "R2": []},
        "denatured": {"R1": [], "R2": []},
    }
    os.makedirs(config["shapemapper_output"], exist_ok=True)
    for rid, sample in tqdm(
        samples.iterrows(), desc="Preparing path ", total=samples.shape[0]
    ):
        title = config["title_template"].format(**dict(sample))
        seq_file = config["sequences"][sample["sequence"]]
        runs[title] = {}
        if os.path.splitext(seq_file)[-1] in [".tsv", ".csv"]:

            seq_filename = os.path.basename(seq_file).split(".")[0]
            fasta_path = os.path.join(
                config["shapemapper_output"], "sequences", seq_filename + ".fasta"
            )
            fasta_from_tsv(seq_file, fasta_path)
            config["sequences"][sample["sequence"]] = fasta_path
            seq_file = fasta_path
        if config["split_seq"]:
            cur_splitted_refs = gen_splitted_ref(
                seq_file,
                os.path.join(config["shapemapper_output"], "sequences"),
            )
            splitted_refs[title] = cur_splitted_refs
            for seq in cur_splitted_refs.keys():
                fastqs = copy.deepcopy(afastq)
                for condtype in fastqs.keys():
                    for readtype in fastqs[condtype].keys():
                        for folder in config["data_folders"]:

                            cur_glob_path = fastq_input_splitted_pattern.format(
                                sequence=seq,
                                input_path=folder,
                                cond=sample[condtype],
                                read=readtype,
                            )
                            cur_glob_path_rev = fastq_input_splitted_pattern.format(
                                sequence=seq + "_rev",
                                input_path=folder,
                                cond=sample[condtype],
                                read=readtype,
                            )

                            # print(f"PATH: {cur_glob_path}")
                            fqs = glob.glob(cur_glob_path) + glob.glob(
                                cur_glob_path_rev
                            )
                            fastqs[condtype][readtype].extend(fqs)

                runs[title][seq] = fastqs
        else:
            fastqs = copy.deepcopy(afastq)
            for condtype in fastqs.keys():
                for readtype in fastqs[condtype].keys():
                    for folder in config["data_folders"]:

                        cur_glob_path = fastq_input_pattern.format(
                            input_path=folder,
                            cond=sample[condtype],
                            read=readtype,
                        )

                        # print(f"PATH: {cur_glob_path}")
                        fqs = glob.glob(cur_glob_path)
                        fastqs[condtype][readtype].extend(fqs)
            runs[title] = fastqs
    return runs, splitted_refs


#    reference: str,
#    output_dir: str,
#    name: str,
#    modified_r1: [str],
#    modified_r2: [str],
#    untreated_r1: [str],
#    untreated_r2: [str],
#    denatured_r1: [str],
#    denatured_r2: [str],
#    cores: int = 8,
#    indiv_norm: bool = True,
#    min_depth: int = 5000,
#    log: str = "log.log",
#    overwrite=True,
#    extra_args: [str] = [],
def launch_shapemapper(config_path, samples_path, interactive=False):
    samples = pd.read_csv(samples_path, sep="\t")
    samples = samples[samples["discard"] != "yes"]
    with open(config_path, "r") as config_file:
        config = YAML.load(config_file)

    runs, splitted_refs = prepare_launch(config, samples)

    check_all_valid = True
    for (title, seqs) in runs.items():
        for seq, fastqs in seqs.items():
            for condtype, strands in fastqs.items():
                if len(strands["R1"]) != len(strands["R2"]):
                    print(f"{title} - {seq} - {condtype} : len(R1) != len(R2)")
                    check_all_valid = False

                for strandname, strand in strands.items():
                    if len(strand) == 0:
                        print(
                            f"{title} - {seq} - {condtype} - {strandname} : no input files"
                        )
                        check_all_valid = False
    if not check_all_valid:
        print("WARNING some conditions have no data")
        if interactive:
            cont = "ask"
            while cont not in ["", "Y", "y", "n", "N", "yes", "no"]:
                cont = input("Would you like to continue [Y/n]")
            if cont == "n" or cont == "N" or cont == "no":
                exit(1)

    for title, seqs in tqdm(
        runs.items(), total=len(runs), desc="Running Shapemapper ", position=0
    ):
        for seq, fastqs in tqdm(
            seqs.items(), total=len(seqs), desc=title, position=1, leave=False
        ):
            run_shapemapper(
                reference=splitted_refs[title][seq],
                output_dir=os.path.join(config["shapemapper_output"], title),
                name=title,
                log=os.path.join(
                    config["shapemapper_output"],
                    title,
                    title + f"_{seq}.log",
                ),
                input_fastqs=fastqs,
            )


def main():
    fire.Fire(launch_shapemapper)


if __name__ == "__main__":
    main()
