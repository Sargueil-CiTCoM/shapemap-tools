import os
import glob
from ruamel.yaml import YAML
import pandas as pd
from . import fasta

yaml = YAML()
base_path = "{path}/{condition}/"
compare_base_path = (
    "{path}/{condition1}/{comp_prefix}{condition1}_{condition2}/"
)
profile_pattern = base_path + "{condition}_{seqid}_profile.txt"
shape_pattern = base_path + "{condition}_{seqid}.shape"
map_pattern = base_path + "{condition}_{seqid}.map"
svg_structure_pattern = base_path + "{condition}_{seqid}{index}.svg"
varna_structure_pattern = base_path + "{condition}_{seqid}.varna"
title_pattern = "{seqid}_{condition}"

aggregate_pattern = base_path + "{condition}_{seqid}_aggregated.tsv"
plot_pattern = base_path + "{condition}_{seqid}_aggregated.svg"
plot_full_pattern = base_path + "{condition}_{seqid}_aggregated_full.svg"
plot_histo_pattern = base_path + "{condition}_{seqid}_histograms.pdf"
plot_histo_svg_pattern = base_path + "{condition}_{seqid}_histograms.svg"
plot_profiles_pattern = base_path + "{condition}_{seqid}_profiles.pdf"
plot_profiles_svg_pattern = base_path + "{condition}_{seqid}_profiles.svg"

# delta_comp_pattern = compare_base_path + "{condition1}{conditions_separator}{condition2}_footprint.svg"
plot_footprint_svg_pattern = compare_base_path + "{condition1}_{condition2}_footprint.svg"
plot_footprint_diff_svg_pattern = compare_base_path + "{condition1}_{condition2}_footprint_diff.svg"
svg_structure_foorprint_pattern = compare_base_path + "{condition1}_{condition2}_footprint_structure.svg"


class ReactivityThreshold:
    INVALID = -0.3
    LOW = 0.40
    MEDIUM = 0.85
    HIGH = 1.0

    COLOR_INVALID = "grey"
    COLOR_NONE = "white"
    COLOR_LOW = "yellow"
    COLOR_MEDIUM = "orange"
    COLOR_HIGH = "red"


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


def get_sequences(path, conditions):
    sequences = set()

    for condition in conditions:
        condlen = len(condition) + 1
        pattern = shape_pattern.format(path=path, condition=condition, seqid="*")

        sequences = sequences.union(
            {
                os.path.splitext(os.path.basename(path))[0][condlen:]
                for path in glob.glob(pattern)
            }
        )
    return sorted(tuple(sequences))


def conditions_from_path(project_path, condition_prefix):
    return sorted(
        [
            os.path.basename(path)
            for path in glob.glob(f"{project_path}/{condition_prefix}*")
        ]
    )


def comparisons_from_path(
    project_path, condition_prefix, comp_prefix="comp_", conditions_separator="__"
):
    return [
        (cond[len(condition_prefix) :] for cond in path.split(conditions_separator))
        for path in glob.glob(f"{project_path}/{condition_prefix}/{comp_prefix}*")
    ]


class Config:
    def __init__(self, config_path):
        with open(config_path, "r") as config_file:
            self.config = yaml.load(config_file)

        sample_path = os.path.join(
            os.path.dirname(config_path), self.config["samples_file"]
        )
        self.samples = pd.read_csv(sample_path, sep="\t")
        self.samples = self.samples[self.samples["discard"] != "yes"]
        self.sequences_id = set()
        # print(self.config["sequences"])
        for seqs_id, seqs_path in self.config["sequences"].items():
            df = pd.read_csv(seqs_path, sep="\t")
            self.sequences_id = self.sequences_id.union(df["name"])

    def parameters(self):
        return list(self.config["parameters"])

    def conditions(self, replicate_col=False):
        ext_params = self.parameters()
        ext_params.append("sequence")
        if replicate_col:
            ext_params.append("replicate")

        res = self.samples[ext_params].drop_duplicates()
        return res

    def path_config(self, reps=["rep0", "rep1", "rep2", "rep3"], input_path=""):
        pc = PathConfig()
        # input_path = (
        #    self.config["shapemapper_output_norm"] if input_path is None else input_path
        # )
        pc.paths = {idx: input_path for idx in reps}

        pc.conditions = {}
        for rid, condition in self.conditions().iterrows():
            cond_dict = condition.to_dict()
            # print("COND", condition.to_dict())
            name = self.config["title_template"].format(replicate="all", **cond_dict)
            pc.conditions[name] = {}
            mcond = self.samples.loc[
                (self.samples.loc[:, cond_dict.keys()] == cond_dict.values()).all(
                    axis=1
                )
            ]
            for rid, subcond in mcond.iterrows():
                pc.conditions[name][subcond["replicate"]] = self.config[
                    "title_template"
                ].format(**subcond)
        pc.sequences = self.sequences_id
        return pc


class PathConfig:
    def __init__(self, path=None):
        if path:
            self._config = self.read_config(path)
            self.paths = self._config["paths"]
            self.conditions = self._config["conditions"]
            self.sequences = self.get_reps_sequences()
        else:
            self.paths = []
            self.conditions = []

    def read_config(self, config_path: str):
        with open(config_path, "r") as fd:
            config = yaml.load(fd)

            for path in config["paths"]:
                os.path.exists(path)

            for condname, condvals in config["conditions"].items():
                for cond_rep_name, cond_rep_path in condvals.items():
                    os.path.exists(
                        os.path.join(config["paths"][cond_rep_name], cond_rep_path)
                    )
            return config

    def get_reps_sequences(self):
        sequences = set()
        for repname, reppath in self.paths.items():
            conds = [
                repcond
                for condname, condreps in self.conditions.items()
                for currepname, repcond in condreps.items()
                if currepname == repname
            ]
            sequences = sequences.union(get_sequences(reppath, conds))
        return sequences

