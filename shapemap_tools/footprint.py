#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import fire


def ratio_sig_test(footprint, ttest_pvalue_thres=0.01, diff_thres=0.2, ratio_thres=0.2):
    footprint = pd.DataFrame(footprint)
    footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
        lambda row: row["delta"]
        if row["pvalue"] <= ttest_pvalue_thres
        and np.abs(row["delta"]) >= diff_thres
        and row["ratio"] >= ratio_thres
        else np.NaN,
        axis=1,
    )
    return footprint


def zfactor_sig_test(footprint, ttest_pvalue_thres=0.01, zfactor_thres=0):
    footprint = pd.DataFrame(footprint)
    footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
        lambda row: row["delta"]
        if row["pvalue"] <= ttest_pvalue_thres and row["zfactor"] > zfactor_thres
        else np.NaN,
        axis=1,
    )
    return footprint


def ttest_only_sig_test(footprint, ttest_pvalue_thres=0.01):
    footprint = pd.DataFrame(footprint)
    footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
        lambda row: row["delta"] if row["pvalue"] <= ttest_pvalue_thres else np.NaN,
        axis=1,
    )
    return footprint


def ratio_zfactor_sig_test(
    footprint,
    ttest_pvalue_thres=0.01,
    diff_thres=0.2,
    ratio_thres=0.2,
    zfactor_thres=0,
):
    footprint = pd.DataFrame(footprint)
    footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
        lambda row: row["delta"]
        if row["pvalue"] <= ttest_pvalue_thres
        and row["zfactor"] > zfactor_thres
        and np.abs(row["delta"]) >= diff_thres
        and row["ratio"] >= ratio_thres
        else np.NaN,
        axis=1,
    )
    return footprint


def footprint_ttest(
    cond1_path,
    cond1_name,
    cond2_path,
    cond2_name,
    deviation_type="stdev",
    zfactor_nsigma=3,
):
    cols = [
        "seqNum",
        "sequence",
        "stat",
        "pvalue",
        "delta",
        "ratio",
        "zfactor",
    ]
    indexes = ["seqNum", "sequence"]
    cond1 = pd.read_csv(cond1_path, sep="\t").set_index(["seqNum", "sequence"])
    cond2 = pd.read_csv(cond2_path, sep="\t").set_index(["seqNum", "sequence"])
    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2}, names=["condition"], axis=1
    )
    res = pd.DataFrame([], columns=cols).set_index(indexes)
    for index, row in footprint.iterrows():

        if row.loc[cond1_name]["desc"] in ("accepted", "reduced") and row.loc[
            cond2_name
        ]["desc"] in ("accepted", "reduced"):

            stat, pvalue = scipy.stats.ttest_ind_from_stats(
                row.loc[cond1_name]["mean"],
                row.loc[cond1_name][deviation_type],
                row.loc[cond1_name]["used_values"],
                row.loc[cond2_name]["mean"],
                row.loc[cond2_name][deviation_type],
                row.loc[cond2_name]["used_values"],
                equal_var=True,
                alternative="two-sided",
            )
            delta = row.loc[cond2_name]["mean"] - row.loc[cond1_name]["mean"]
            ratio = np.abs(delta) / (
                row.loc[cond2_name]["mean"] + row.loc[cond1_name]["mean"]
            )
            try:
                zfactor = (
                    1
                    - (
                        zfactor_nsigma
                        * (
                            row.loc[cond2_name][deviation_type]
                            + row.loc[cond1_name][deviation_type]
                        )
                    )
                    / delta
                )
            except ZeroDivisionError: 
                print("Z-factor computation failed : DivisionByZero")

            curres = pd.DataFrame(
                [[index[0], index[1], stat, pvalue, delta, ratio, zfactor]],
                columns=cols,
            ).set_index(indexes)
            res = pd.concat([res, curres])
        else:
            curres = pd.DataFrame(
                [[index[0], index[1], np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]],
                columns=cols,
            ).set_index(indexes)
            res = pd.concat([res, curres])

    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2, "ttest": res},
        names=["condition"],
        axis=1,
    )

    return footprint


def plot_reactivity(
    footprint: pd.DataFrame,
    title="Footprint",
    output="fig.svg",
    format="svg",
    deviation_type="stdev",
    cond1_name="Condition1",
    cond2_name="Condition2",
):

    replicates = footprint.drop("ttest", axis=1)

    means = replicates.xs("mean", level=1, axis=1).reset_index()

    means = means.replace(-10, np.NaN)

    unidmeans = means.drop(["seqNum", "sequence"], axis=1)

    ax = unidmeans[unidmeans.columns[0]].plot(
        kind="bar",
        width=1.0,
        align="center",
        figsize=(len(footprint) / 3, 4),
        # alpha=0.5,
        color="green",
    )
    unidmeans[unidmeans.columns[1]].plot(
        kind="bar",
        width=1.0,
        align="center",
        figsize=(len(footprint) / 3, 4),
        alpha=0.5,
        color="lime",
        ax=ax,
    )

    ttest = footprint["ttest"]

    ttest = ttest.reset_index()

    ttest["xlabel"] = ttest["seqNum"].astype(str) + "\n" + ttest["sequence"].astype(str)
    ttest["color"] = "white"
    ttest.loc[ttest["significant_delta"] > 0, "color"] = "yellow"
    ttest.loc[ttest["significant_delta"] < 0, "color"] = "orange"

    first_idx = ttest.first_valid_index()
    last_idx = ttest.last_valid_index()

    for index, row in ttest.loc[first_idx:last_idx].iterrows():
        idx = index
        if row["color"] != "white":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color=row["color"],
            )

    ax.axhline(y=0.4, color="orange", linestyle="-", label="Medium Reactivity")
    ax.axhline(y=0.7, color="red", linestyle="-", label="High reactivity")

    plt.margins(0)
    legend_elements = [
        Patch(facecolor="green", label=f"{cond1_name}"),
        Patch(facecolor="lime", label=f"{cond2_name}", alpha=0.5),
        Patch(facecolor="yellow", label=f"{cond1_name} sign. Higher", alpha=0.3),
        Patch(facecolor="orange", label=f"{cond2_name} sign. Higher", alpha=0.3),
        Line2D([0], [0], color="red", label="High reactivity threshold"),
        Line2D([0], [0], color="orange", label="Medium reactivity threshold"),
    ]
    plt.legend(handles=legend_elements, loc="upper left")
    plt.title(title, loc="left")
    plt.tight_layout()
    plt.savefig(output, format=format)


# def footprint_ttest(
#        cond1_path,
#        cond1_name,
#        cond2_path,
#        cond2_name,
#        sigtest=ratio_sig_test,
#        deviation_type="mad",
#        ttest_pvalue_thres=0.01,
#        diff_thres=0.2,
#        ratio_thres=0.2,
#        zfactor_nsigma=3,
#        ):


def footprint_main(
    cond1: str,
    cond2: str,
    cond1_name: str = None,
    cond2_name: str = None,
    deviation_type: str = "stdev",
    ttest_pvalue_thres=0.01,
    diff_thres=0.2,
    ratio_thres=0.2,
    output: str = None,
    plot: str = None,
    plot_title="Footprint",
    plot_format="svg",
):
    cond1_name = cond1_name if cond1_name is not None else cond1
    cond2_name = cond2_name if cond2_name is not None else cond2
    assert deviation_type in ["stdev", "sem", "mad"]

    footprint = footprint_ttest(
        cond1,
        cond1_name,
        cond2,
        cond2_name,
        deviation_type=deviation_type,
    )

    footprint = ratio_sig_test(
        footprint,
        ttest_pvalue_thres=ttest_pvalue_thres,
        diff_thres=diff_thres,
        ratio_thres=ratio_thres,
    )

    if output is not None:
        footprint.to_csv(output, sep="\t")

    if plot is not None:
        plot_reactivity(
            footprint,
            plot_title,
            plot,
            plot_format,
            deviation_type=deviation_type,
            cond1_name=cond1_name,
            cond2_name=cond2_name,
        )


def main():
    return fire.Fire(footprint_main)


if __name__ == "__main__":
    main()
