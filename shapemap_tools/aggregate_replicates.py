#!/usr/bin/env python
# coding: utf-8

# In[1]:

import glob
import pandas as pd
import numpy as np
import os
import copy
from matplotlib import cm
import matplotlib.pyplot as plt
import gc
import fire
from . import utils

# pd.set_option("display.max_rows", None)  # or 1000


def plot_aggregate(
    aggregated: pd.DataFrame,
    fulloutput="fig.full.svg",
    output="fig.svg",
    title: str = "Aggregated reactivity",
    format="svg",
):
    aggregated = copy.deepcopy(aggregated)
    aggregated["xlabel"] = (
        aggregated.index.get_level_values("seqNum").astype(str)
        + "\n"
        + aggregated.index.get_level_values("sequence").astype(str)
    )
    replicates = aggregated.loc[
        :,
        aggregated.columns.drop(["mean", "stdev", "sem", "mad"]),
    ].replace(-10, np.NaN)

    aggregated = aggregated.sort_index()
    meanstdev = aggregated.loc[:, ["xlabel", "mean", "stdev"]].replace(-10, np.NaN)

    reactivities_cols = [col for col in aggregated.columns if col[1] == "reactivity"]
    reactivities_cols += [("xlabel", "")]

    ax = replicates[reactivities_cols].plot(
        x="xlabel",
        kind="bar",
        width=0.7,
        stacked=False,
        figsize=(len(aggregated) / 3.5, 4),
        align="center",
        xticks=np.arange(0, len(aggregated) + 1, 1),
    )
    ax = meanstdev[["mean", "xlabel"]].plot(
        x="xlabel",
        y="mean",
        drawstyle="steps-mid",
        ax=ax,
        colormap=cm.cubehelix,
        linewidth=0.5,
    )
    ax.set_xlabel("Position")
    ax.set_ylabel("Reactivity")
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.margins(0)
    plt.title(title, loc="left")
    plt.legend(loc="upper left")
    ax.errorbar(
        meanstdev.index.get_level_values("seqNum") - 1,
        meanstdev["mean"],
        yerr=meanstdev["stdev"],
        fmt="",
        color="k",
        ls="none",
        capsize=4,
    linewidth=0.5,
    )
    try:
        plt.tight_layout()
        # plt.show()
        plt.savefig(fulloutput, format=format)
    except ValueError as e:
        print(f"Unable to save fullplot: {e}")
        open(fulloutput, "a").close()

    fig, ax = plt.subplots(
        nrows=1, ncols=1, sharex=True, figsize=(len(aggregated) / 4, 4)
    )

    aggregated["color"] = utils.ReactivityThreshold.COLOR_NONE
    aggregated.loc[
        (aggregated["mean"] > utils.ReactivityThreshold.HIGH), "color"
    ] = utils.ReactivityThreshold.COLOR_HIGH
    aggregated.loc[
        (
            (aggregated["mean"] <= utils.ReactivityThreshold.HIGH)
            & (aggregated["mean"] > utils.ReactivityThreshold.MEDIUM)
        ),
        "color",
    ] = utils.ReactivityThreshold.COLOR_MEDIUM
    aggregated.loc[
        (
            (aggregated["mean"] <= utils.ReactivityThreshold.MEDIUM)
            & (aggregated["mean"] > utils.ReactivityThreshold.LOW)
        ),
        "color",
    ] = utils.ReactivityThreshold.COLOR_LOW
    aggregated.loc[
        (aggregated["mean"] < utils.ReactivityThreshold.INVALID), "color"
    ] = utils.ReactivityThreshold.COLOR_INVALID

    aggregated.loc[(aggregated["mean"] == -10), "stdev"] = np.NaN
    aggregated.loc[(aggregated["mean"] == -10), "mean"] = np.NaN
    aggregated["xlabel_rot"] = (
        aggregated.index.get_level_values("seqNum").astype(str)
        + " - "
        + aggregated.index.get_level_values("sequence").astype(str)
    )

    aggregated.plot(
        ax=ax,
        # x="xlabel_rot",
        rot=70,
        y="mean",
        kind="bar",
        width=1,
        color=aggregated["color"],
        yerr="stdev",
        stacked=False,
        capsize=3,
    )
    ax = meanstdev[["mean", "xlabel"]].plot(
        x="xlabel",
        y="mean",
        drawstyle="steps-mid",
        ax=ax,
        colormap=cm.cubehelix,
        linewidth=0.5,
    )

    ax.set_xlabel("Position")
    ax.set_ylabel("Reactivity")
    plt.margins(0)
    plt.title(title, loc="left")
    plt.legend(loc="upper left")
    try:
        plt.tight_layout()
        # plt.show()
        plt.savefig(output, format=format)
    except ValueError as e:
        print(f"Unable to save plot: {e}")
        open(output, "a").close()

def aggregate_replicates(
    output_path: str = "output",
    input_dirs: [str] = [],
    path_config_path: str = None,
    config_path: str = None,
):
    if path_config_path:
        config = utils.PathConfig(path_config_path)
    else:
        base_config = utils.Config(config_path)
        config = base_config.path_config()

    sequences = config.sequences
    assert len(sequences) > 0
    #print(sequences)
    profiles = {}

    for condname, reps in config.conditions.items():
        print(f"Aggregating {condname}")
        profiles[condname] = {}
        for seqid in sequences:
            repsdf = {}
            for rep_id, rep_path in reps.items():
                try:
                    print(rep_path)
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
                profiles[condname][seqid]["mad"] = profiles[condname][seqid][conds].mad(
                    axis=1
                )
            else:
                profiles[condname][seqid] = None

    os.makedirs(output_path, exist_ok=True)
    for condname, seqs_profiles in profiles.items():
        os.makedirs(
            utils.base_path.format(path=output_path, condition=condname), exist_ok=True
        )
        for seqid, profile in seqs_profiles.items():
            if profile is not None:
                print(profile)
                profile.to_csv(
                    utils.aggregate_pattern.format(
                        path=output_path, condition=condname, seqid=seqid
                    ),
                    sep="\t",
                    
                )
                write_shape(output_path, condname, seqid, profile)
                write_map(output_path, condname, seqid, profile)
    #plot_all(output_path, profiles)


def write_shape(output_path, condition, seqid, profile):
    shprofile = profile.reset_index(level="sequence")[["mean"]]
    idxmin = shprofile.index.min()
    firstrows = pd.DataFrame({"mean": np.full(idxmin - 1, -10)}, index=range(1, idxmin))
    firstrows.index.names = ["seqNum"]
    shprofile = pd.concat([firstrows, shprofile])
    shprofile.to_csv(
        utils.shape_pattern.format(path=output_path, condition=condition, seqid=seqid),
        sep="\t",
        float_format="%.4f",
        header=False,
    )
    return shprofile


def write_map(output_path, condition, seqid, profile):
    mapprofile = profile.reset_index(level="sequence")[["mean", "stdev", "sequence"]]
    idxmin = mapprofile.index.min()
    firstrows = pd.DataFrame(
        {
            "mean": np.full(idxmin - 1, -10),
            "stdev": np.zeros(idxmin - 1),
            "sequence": np.full(idxmin - 1, "N"),
        },
        index=range(1, idxmin),
    )
    firstrows.index.names = ["seqNum"]
    mapprofile = pd.concat([firstrows, mapprofile])
    mapprofile.to_csv(
        utils.map_pattern.format(path=output_path, condition=condition, seqid=seqid),
        sep="\t",
        float_format="%.4f",
        header=False,
    )
    return mapprofile


def plot_all(output_path, profiles):
    for cond, seqs_profiles in profiles.items():
        print(f"Plotting {cond}")
        for seqid, curdf in seqs_profiles.items():
            print(f"Plotting {seqid}")
            if curdf is not None:
                plot_aggregate(
                    curdf,
                    fulloutput=utils.plot_full_pattern.format(
                        path=output_path, condition=cond, seqid=seqid
                    ),
                    output=utils.plot_pattern.format(
                        path=output_path, condition=cond, seqid=seqid
                    ),
                    title=f"{cond} - {seqid} Reactivity",
                )
                gc.collect()
            # curdf.to_csv(f"{cond}/{cond}_{aptaid}_aggregated.tsv",sep="\t")
            # curdf = curdf.reset_index(drop=False)
            # curdf["mean"] = curdf["mean"].fillna(-999)
            # curdf["stdev"] = curdf["stdev"].fillna(0)
            # curdf[["seqNum", "mean", "stdev", "sequence"]]
            # .to_csv(f"{cond}/{cond}_{aptaid}_aggregated.map",sep="\t",
            # index=False, header=None)
            # curdf[["seqNum", "mean"]]
            # .to_csv(f"{cond}/{cond}_{aptaid}_aggregated.shape",sep="\t",
            # header=None,index=False)


def main():
    fire.Fire(aggregate_replicates)


if __name__ == "__main__":
    main()
