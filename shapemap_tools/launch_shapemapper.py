import fire
import subprocess as sp
import pandas as pd
import ruamel.yaml as yaml
import glob
import copy
from . import fasta, utils
import os
from tqdm import tqdm
import string
import random

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
    shapemapper_path: str = "shapemapper2",
    seq="",
    verbose=False,
):
    random_tmp = "".join(random.choice(string.ascii_lowercase) for i in range(8))
    args = [
        "--name",
        name,
        "--nproc",
        str(cores),
        "--min-depth",
        str(min_depth),
        # "--amplicon",
        "--target",
        reference,
        "--log",
        log,
        "--temp",
        f"temp/{name}_{seq}_{random_tmp}",
        "--out",
        output_dir,
    ]
    if indiv_norm:
        args.append("--indiv-norm")
    if overwrite:
        args.append("--overwrite")
    if verbose:
        args.append("--verbose")

    form_inputs = []
    for condtype in input_fastqs.keys():
        form_inputs.append(f"--{condtype}")
        for readtype, fastqs in input_fastqs[condtype].items():
            form_inputs.append(f"--{readtype}")
            form_inputs.extend(fastqs)

    cmd = [shapemapper_path] + form_inputs + args
    # print(" ".join(cmd))
    try:
        sp.run(cmd, capture_output=True, text=True)
    except Exception as e:
        print(f"Error launching shapemapper for {name} {e}")
        print(" ".join(cmd))


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


def get_config(config, path, default):
    cur_conf = config
    for p in path:
        if p in cur_conf:
            cur_conf = cur_conf[p]
        else:
            return default
    return cur_conf


def prepare_launch(config, samples, dnerase=False):
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
            seq_conf = get_config(config, ["sequence_config"], {})
            fasta_from_tsv(
                seq_file,
                fasta_path,
                name_col=get_config(seq_conf, "name_col", "name"),
                seq_col=get_config(seq_conf, "seq_col", "sequence"),
                prefix_cols=get_config(seq_conf, "prefix_cols", []),
                suffix_cols=get_config(seq_conf, "suffix_cols", []),
            )

            config["sequences"][sample["sequence"]] = fasta_path
            seq_file = fasta_path
        if config["split_seq"]:
            cur_splitted_refs = gen_splitted_ref(
                seq_file,
                os.path.join(config["shapemapper_output"], "sequences"),
            )
            splitted_refs[title] = cur_splitted_refs
            for seq in cur_splitted_refs.keys():
                if not dnerase or (
                    dnerase
                    and not os.path.exists(
                        utils.profile_pattern.format(
                            path=config["shapemapper_output"],
                            seqid=seq,
                            condition=title,
                        )
                    )
                ):
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
                        # print(fqs)
                        fastqs[condtype][readtype].extend(fqs)
            runs[title] = {"all": fastqs}
            splitted_refs[title] = {"all": seq_file}
    # print(runs)
    # print(splitted_refs)
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
def launch_shapemapper(
    config_path,
    samples_path,
    interactive=False,
    shapemapper_path=None,
    dnerase=False,
    verbose=False,
    nthreads=8,
):
    samples = pd.read_csv(samples_path, sep="\t")
    samples = samples[samples["discard"] != "yes"]
    with open(config_path, "r") as config_file:
        config = YAML.load(config_file)

    if "shapemapper_path" in config and shapemapper_path is None:
        shapemapper_path = config["shapemapper_path"]

    assert shapemapper_path is not None
    runs, splitted_refs = prepare_launch(config, samples, dnerase=dnerase)

    check_all_valid = True
    for title, seqs in runs.items():
        # print(title)
        # print(seqs)

        for seq, fastqs in seqs.items():
            # print(seq)
            # print(fastqs)
            for condtype, strands in fastqs.items():
                # print(condtype)
                # print(strands)
                if len(strands["R1"]) != len(strands["R2"]):
                    print(
                        f"{title} - {seq} - {condtype} : len(R1) {len(strands['R1'])}!= len(R2) {len(strands['R2'])}"
                    )
                    check_all_valid = False

                for strandname, strand in strands.items():
                    # print(strandname)
                    # print(strand)
                    # print(len(strand))
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
            if not dnerase or (
                dnerase
                and not os.path.exists(
                    utils.profile_pattern.format(
                        path=config["shapemapper_output"], seqid=seq, condition=title
                    )
                )
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
                    shapemapper_path=shapemapper_path,
                    seq=seq,
                    verbose=verbose,
                    cores=nthreads,
                )


def main():
    fire.Fire(launch_shapemapper)


if __name__ == "__main__":
    main()
