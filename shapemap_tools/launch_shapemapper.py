import fire
import subprocess as sp
import pandas as pd
import ruamel.yaml as yaml
import glob
import copy
from . import fasta
import os

YAML = yaml.YAML()

fastq_input_pattern = "{input_path}/{cond}/*{read}*.fastq*"
fastq_input_splitted_pattern = "{input_path}/{cond}/*{read}*_{sequence}.fastq*"


def run_shapemapper(
    reference: str,
    output_dir: str,
    name: str,
    modified_r1: [str],
    modified_r2: [str],
    untreated_r1: [str],
    untreated_r2: [str],
    denatured_r1: [str],
    denatured_r2: [str],
    cores: int = 8,
    indiv_norm: bool = True,
    min_depth: int = 5000,
    log: str = "log.log",
    overwrite=True,
    extra_args: [str] = [],
):
    cmd = [
        "shapemapper",
        "--name",
        "name",
        "--nproc",
        str(cores),
        "--min-depth",
        str(min_depth),
        "--target",
        reference,
        "--out",
        output_dir,
        "--log",
        log,
        "--modified",
        "--R1",
        " ".join(modified_r1),
        "--R2",
        " ".join(modified_r2),
        "--untreated",
        "--R1",
        " ".join(untreated_r1),
        "--R2",
        " ".join(untreated_r2),
        "--denatured",
        "--R1",
        " ".join(denatured_r1),
        "--R2",
        " ".join(denatured_r2),
    ]
    if indiv_norm:
        cmd.append("--indiv-norm")
    if overwrite:
        cmd.append("--overwrite")

    sp.run(cmd)


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
def launch_shapemapper(config_path, samples_path):
    samples = pd.read_csv(samples_path, sep="\t")
    config = YAML.load(config_path)
    afastq = {
        "modified": {"R1": [], "R2": []},
        "untreated": {"R1": [], "R2": []},
        "denatured": {"R1": [], "R2": []},
        }
    os.makedirs(config["shapemapper_output"], exist_ok=True)
    for rid, sample in samples.rows():
        if config["split_seq"]:
            splitted_refs = gen_splitted_ref(
                sample["sequence"],
                os.path.join(config["shapemapper_output"], "sequences"),
            )
            for seq in splitted_refs.keys():
                fastqs = copy.deepcopy(afastq)
                for condtype in fastqs.keys():
                    for readtype in fastqs[condtype].keys():
                        for folder in config["data_folders"]:
                            fqs = glob.glob(
                                fastq_input_splitted_pattern.format(
                                    sequence=seq,
                                    input_path=folder,
                                    cond=sample[condtype],
                                    read=readtype,
                                )
                            )
                            fastqs[condtype][readtype].extend(fqs)
                run_shapemapper(
                    reference=splitted_refs[sample["sequence"]],
                    output_dir=config["shapemapper_output"],
                    name=config["title_template"].format(**dict(sample))
                    + f"_{sample['sequence']}",
                    modified_r1=fastqs["modified"]["R1"],
                    modified_r2=fastqs["modified"]["R2"],
                    untreated_r1=fastqs["untreated"]["R1"],
                    untreated_r2=fastqs["untreated"]["R2"],
                    denatured_r1=fastqs["denatured"]["R1"],
                    denatured_r2=fastqs["denatured"]["R2"],
                )


def main():
    fire.Fire(launch_shapemapper)


if __name__ == "__main__":
    main()
