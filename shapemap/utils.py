import os
import glob
from ruamel.yaml import YAML

yaml = YAML()
base_path = "{path}/{condition}/"
profile_pattern = base_path + "{condition}_{seqid}_profile.txt"
shape_pattern = base_path + "{condition}_{seqid}.shape"
map_pattern = base_path + "{condition}_{seqid}.map"
title_pattern = "{seqid} {condition}"
aggregate_pattern = base_path + "{condition}_{seqid}_aggregated.tsv"
plot_pattern = base_path + "{condition}_{seqid}_aggregated.svg"
plot_full_pattern = base_path + "{condition}_{seqid}_aggregated_full.svg"
plot_histo_pattern = base_path + "{condition}_{seqid}_histograms.pdf"
plot_profiles_pattern = base_path + "{condition}_{seqid}_profiles.pdf"


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
    return sequences


class Config:
    def __init__(self, path):
        self._config = self.read_config(path)
        self.paths = self._config["paths"]
        self.conditions = self._config["conditions"]

        self.sequences = self.get_sequences()

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
