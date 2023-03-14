#!/usr/bin/env python
# coding: utf-8

# In[1]:

# import glob
import pandas as pd
import numpy as np
import os
import fire
from . import utils
from tqdm import tqdm

# pd.set_option("display.max_rows", None)  # or 1000

import warnings


def aggregate_replicates(
    output_path: str = "output",
    input_path: str = "intput",
    path_config_path: str = None,
    config_path: str = None,
):
    if path_config_path:
        config = utils.PathConfig(path_config_path)
    else:
        base_config = utils.Config(config_path)
        # # FIXME find info about replicate in config file
        config = base_config.path_config(input_path=input_path)

    sequences = config.sequences
    assert len(sequences) > 0
    # print(sequences)
    profiles = {}

    for condname, reps in tqdm(
        config.conditions.items(), total=len(config.conditions.items()), position=0
    ):
        profiles[condname] = {}
        for seqid in tqdm(
            sequences,
            total=len(sequences),
            position=1,
            leave=0,
            desc=f"Aggregating {condname}",
        ):
            repsdf = {}
            for rep_id, rep_path in reps.items():
                try:
                    # print(rep_path)
                    curdf = pd.read_csv(
                        utils.map_pattern.format(
                            path=config.paths[rep_id],
                            condition=rep_path,
                            seqid=seqid,
                        ),
                        sep="\t",
                        names=["seqNum", "reactivity", "sterr", "sequence"],
                    )
                    curdf = curdf.set_index(["seqNum", "sequence"])
                    curdf = curdf.replace(-999, np.NaN)
                    repsdf[rep_id] = curdf
                    # profiles[condname][seqid][rep_id] = curdf["reactivity"]

                except FileNotFoundError as fnfe:
                    # print(lol)
                    print(fnfe)
            if len(repsdf) > 0:
                profiles[condname][seqid] = pd.concat(repsdf, axis=1)
                conds = [
                    col
                    for col in profiles[condname][seqid].columns
                    if col[1] == "reactivity"
                ]
                profiles[condname][seqid]["mean"] = profiles[condname][seqid][
                    conds
                ].mean(axis=1)
                profiles[condname][seqid]["stdev"] = profiles[condname][seqid][
                    conds
                ].std(axis=1, ddof=1)
                profiles[condname][seqid]["median"] = profiles[condname][seqid][
                    conds
                ].median(axis=1)
                profiles[condname][seqid]["sem"] = profiles[condname][seqid][conds].sem(
                    axis=1, ddof=1
                )

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    profiles[condname][seqid]["mad"] = profiles[condname][seqid][
                        conds
                    ].mad(axis=1)
            else:
                profiles[condname][seqid] = None

    os.makedirs(output_path, exist_ok=True)
    for condname, seqs_profiles in tqdm(
        profiles.items(), total=len(profiles.items()), desc="Writing files"
    ):
        os.makedirs(
            utils.base_path.format(path=output_path, condition=condname), exist_ok=True
        )
        for seqid, profile in seqs_profiles.items():
            if profile is not None:
                # print(profile)
                profile.to_csv(
                    utils.aggregate_pattern.format(
                        path=output_path, condition=condname, seqid=seqid
                    ),
                    sep="\t",
                )
                profile["mean"] = profile["mean"].fillna(-999)
                write_shape(output_path, condname, seqid, profile)
                write_map(output_path, condname, seqid, profile)
    # plot_all(output_path, profiles)


def write_shape(output_path, condition, seqid, profile):
    shprofile = profile.reset_index(level="sequence")[["mean"]]
    shprofile.to_csv(
        utils.shape_pattern.format(path=output_path, condition=condition, seqid=seqid),
        sep="\t",
        float_format="%.4f",
        header=False,
    )
    return shprofile


def write_map(output_path, condition, seqid, profile):
    mapprofile = profile.reset_index(level="sequence")[["mean", "stdev", "sequence"]]

    mapprofile.to_csv(
        utils.map_pattern.format(path=output_path, condition=condition, seqid=seqid),
        sep="\t",
        float_format="%.4f",
        header=False,
    )
    return mapprofile


def main():
    fire.Fire(aggregate_replicates)


if __name__ == "__main__":
    main()
