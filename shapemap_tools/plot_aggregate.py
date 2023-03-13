#!/usr/bin/env python
# coding: utf-8

# In[1]:

import glob
import pandas as pd
import numpy as np
import os
import copy
import multiprocessing as mp
from matplotlib import cm
import matplotlib.pyplot as plt
import gc
import fire
from . import utils
from tqdm import tqdm

import warnings

# pd.set_option("display.max_rows", None)  # or 1000


def plot_aggregate_wrapper(args):
    plot_aggregate(**args)


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

    meanstdev = (
        aggregated.loc[:, ["xlabel", "mean", "stdev"]]
        .replace(-10, np.NaN)
        .droplevel(1, axis=1)
    )
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
    # print(meanstdev)
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

    aggregated = aggregated.droplevel(1, axis=1)
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


def plot_aggregates(path: str, dnerase: bool = False, nthreads=1):
    files = glob.glob(f"{path}/**/*_aggregated.tsv") + glob.glob(
        f"{path}/*_aggregated.tsv"
    )

    tasks = []
    for file in tqdm(files, total=len(files)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                profile = pd.read_csv(file, sep="\t", header=[0, 1], index_col=[0, 1])
                out = os.path.splitext(file)[0] + ".svg"
                fullout = os.path.splitext(file)[0] + "_full.svg"
                if (
                    not dnerase
                    or not os.path.exists(out)
                    or not os.path.exists(fullout)
                ):
                    if nthreads <= 1:
                        plot_aggregate(
                            profile,
                            output=out,
                            fulloutput=fullout,
                            title=f"{os.path.splitext(os.path.basename(file))[0]} Reactivity",
                        )
                    else:
                        tasks += [
                            {
                                "aggregated": profile,
                                "output": out,
                                "fulloutput": fullout,
                                "title": f"{os.path.splitext(os.path.basename(file))[0]} Reactivity",
                            }
                        ]
            except Exception as e:
                print(f"Error Plotting {file}")
                print(e)
            gc.collect()
    if nthreads > 1:
        with mp.Pool(nthreads) as pool:
            tqdm(pool.imap_unordered(plot_aggregate_wrapper, tasks), total=len(tasks))


def main():
    fire.Fire(plot_aggregates)


if __name__ == "__main__":
    main()
