import fire
import subprocess as sp
import pandas as pd
import ruamel.yaml as yaml
import glob
import copy
from . import fasta
import os

shapemapper_path = "/data/fxlyonnet/shapemapper-2.1.5/shapemapper"

YAML = yaml.YAML()

fastq_input_pattern = "{input_path}/{cond}/*{read}*.fastq*"
fastq_input_splitted_pattern = "{input_path}/{cond}/*{read}*_{sequence}.fastq*"


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

    print(" ".join(cmd))
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


def prepare_launch(config, samples):
    runs = {}
    splitted_refs = {}
    afastq = {
        "modified": {"R1": [], "R2": []},
        "untreated": {"R1": [], "R2": []},
        "denatured": {"R1": [], "R2": []},
    }
    os.makedirs(config["shapemapper_output"], exist_ok=True)
    for rid, sample in samples.iterrows():
        title = config["title_template"].format(**dict(sample))
        runs[title] = {}
        if config["split_seq"]:
            cur_splitted_refs = gen_splitted_ref(
                config["sequences"][sample["sequence"]],
                os.path.join(config["shapemapper_output"], "sequences"),
            )
            splitted_refs[title] = cur_splitted_refs
            for seq in cur_splitted_refs.keys():
                #print(f"Preparing path for {title} - {seq}")
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
def launch_shapemapper(config_path, samples_path):
    samples = pd.read_csv(samples_path, sep="\t")
    samples = samples[samples["discard"] != "yes"]
    with open(config_path, "r") as config_file:
        config = YAML.load(config_file)

    runs, splitted_refs = prepare_launch(config, samples)

    check_all_valid = True
    for (title, seqs) in runs.items():
        for seq, fastqs in seqs.items():
            for condtype, strands in fastqs.items():
                for strandname, strand in strands.items():
                    if len(strand) == 0:
                        print(f"{title} - {seq} - {condtype} - {strandname} : no input files")
                        check_all_valid = False
    if not check_all_valid:
        print("WARNING some conditions have no data") 

    for title, seqs in runs.items():
        for seq, fastqs in seqs.items():
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
