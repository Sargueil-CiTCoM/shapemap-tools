#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import fire
from tqdm import tqdm
import os

# def ratio_sig_test(footprint, ttest_pvalue_thres=0.01, diff_thres=0.2, ratio_thres=0.2):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"]
#         if row["pvalue"] <= ttest_pvalue_thres
#         and np.abs(row["delta"]) >= diff_thres
#         and row["ratio"] >= ratio_thres
#         else np.NaN,
#         axis=1,
#     )
#     return footprint


# def zfactor_sig_test(footprint, ttest_pvalue_thres=0.01, zfactor_thres=0):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"]
#         if row["pvalue"] <= ttest_pvalue_thres and row["zfactor"] > zfactor_thres
#         else np.NaN,
#         axis=1,
#     )
#     return footprint


# def ttest_only_sig_test(footprint, ttest_pvalue_thres=0.01):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"] if row["pvalue"] <= ttest_pvalue_thres else np.NaN,
#         axis=1,
#     )
#     return footprint


# def ratio_zfactor_sig_test(
#     footprint,
#     ttest_pvalue_thres=0.01,
#     diff_thres=0.2,
#     ratio_thres=0.2,
#     zfactor_thres=0,
# ):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"]
#         if row["pvalue"] <= ttest_pvalue_thres
#         and row["zfactor"] > zfactor_thres
#         and np.abs(row["delta"]) >= diff_thres
#         and row["ratio"] >= ratio_thres
#         else np.NaN,
#         axis=1,
#     )
#     return footprint


def footprint_ttest(
    cond1_path,
    cond1_name,
    cond2_path,
    cond2_name,
    deviation_type="stdev",
    # zfactor_nsigma=3,
    ttest_pvalue_thres=0.05,
    diff_thres=0.2,
    ratio_thres=0.2,
):
    """
    Compares two reactivity profiles (two 'conditions')

    Args:
        cond1_path: aggreact_tsv file of condition 1
        cond1_name: name of condition 1
        cond2_path: aggreact_tsv file of condition 2
        cond1_name: name of condition 2
        
    Returns:
        footprint (DataFrame): results prepared for plotting by footprint.plot_reactivity()
        footprint_csv (DataFrame): results in table format (to be written to tsv file)
    """

    cols = [
        "seqNum",
        "sequence",
        "pvalue",
        "difference",
        "ratio",
        "significant_higher",
        "significant_lower"
    ]
    indexes = ["seqNum", "sequence"]

    df = pd.read_csv(cond1_path, sep="\t", header=None, nrows=3)
    row1 = df.iloc[0].tolist()
    row2 = df.iloc[1].tolist()
    row3 = df.iloc[2].tolist()
    new_columns = row3[:2] + row2[2:-4] +row1[-4:]
    cond1 = pd.read_csv(cond1_path, sep="\t", header=None, skiprows=3)
    cond2 = pd.read_csv(cond2_path, sep="\t", header=None, skiprows=3)
    cond1.columns = new_columns
    cond2.columns = new_columns
    cond1 = cond1.set_index(["seqNum", "sequence"])
    cond2 = cond2.set_index(["seqNum", "sequence"])
    
    #joins rows of the two aggreact tables with same seqNum (and nucleotide)
    # (WARNING: if cond2 starts at lower seqNum, these positions are at the end!)
    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2}, names=["condition"], axis=1
    )
    res = list()
    for index, row in footprint.iterrows():
        if np.isnan(row.loc[cond1_name]["mean"]) == False and np.isnan(row.loc[
            cond2_name]["mean"]) == False:
            # get reactivities of the replicates
            rs1 = np.array(list(row.loc[cond1_name]['reactivity']))
            rs2 = np.array(list(row.loc[cond2_name]['reactivity']))
            # compare by ttest only if at least one condition has more than
            # 1 replicates
            if len(rs1[~np.isnan(rs1)])>1 or len(rs2[~np.isnan(rs2)])>1:
                _, pvalue = scipy.stats.ttest_ind(
                    rs1, rs2,
                    equal_var=True,
                    alternative="two-sided",nan_policy='omit')
            else:
                pvalue = 0 # set to 0 in order to leave decision to other criteria

            difference = row.loc[cond2_name]["mean"] - row.loc[cond1_name]["mean"]
            mean_sum = row.loc[cond2_name]["mean"] + row.loc[cond1_name]["mean"]
            if mean_sum!=0:
                # check again: what is the exact idea of ratio??
                ratio = np.abs(difference) / (mean_sum)
            else: ratio = 0

            significant_higher = 'NO'
            significant_lower = 'NO'
            
            if pvalue < ttest_pvalue_thres and ratio > ratio_thres:
                if difference > diff_thres:
                   significant_higher = 'YES'
                elif difference < -diff_thres:
                   significant_lower = 'YES'

            curres = pd.DataFrame(
                [[index[0], index[1], pvalue, difference, ratio, significant_higher, significant_lower]],
                columns=cols,).set_index(indexes)
            res.append(curres)
        else:
            curres = pd.DataFrame(
                [[index[0], index[1], np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]],
                columns=cols,).set_index(indexes)
            res.append(curres)

    #res = pd.DataFrame([], columns=cols).set_index(indexes)
    res = pd.concat(res)

    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2, "ttest": res},
        names=["condition"],
        axis=1,)

    footprint_csv = pd.concat(
        {cond1_name: cond1["mean"], cond2_name: cond2["mean"], "analysis": res},
        axis=1)
    footprint_csv = footprint_csv.sort_values(by=["seqNum"])

    return footprint, footprint_csv


def plot_reactivity(
    footprint: pd.DataFrame,
    title="Footprint",
    dif_title="Footprint",
    format="svg",
    output="fig.svg",
    diff_output="diff.svg",
    cond1_name="Condition1",
    cond2_name="Condition2",
):

    replicates = footprint.drop("ttest", axis=1)
    means = replicates.xs("mean", level=1, axis=1).reset_index()
    means = means.sort_values(by=["seqNum"]).reset_index().drop(["index"], axis=1)
    means["xlabel"] = (means["seqNum"].astype(str) + "\n" + means["sequence"].astype(str))
    means = means.replace(-10, -0.1)
    unidmeans = means.drop(["seqNum", "sequence"], axis=1)

    stdev = replicates.xs("stdev", level=1, axis=1).reset_index()
    stdev = stdev.sort_values(by=["seqNum"]).reset_index().drop(["index"], axis=1)
    stdev["xlabel"] = (stdev["seqNum"].astype(str) + "\n" + stdev["sequence"].astype(str))
    unidstdev = stdev.drop(["seqNum", "sequence"], axis=1)

    significant = footprint[[('ttest', 'significant_higher'),('ttest', 'significant_lower')]].reset_index()
    significant = significant.sort_values(by=[('seqNum', '')]).reset_index().drop([('index','')], axis=1)
    significant["xlabel"] = (significant[('seqNum', '')].astype(str) + "\n" + significant[('sequence','')].astype(str))
    unidsignificant = significant.drop([('seqNum', ''), ('sequence','')], axis=1)
   
    difference = footprint[('ttest', 'difference')].reset_index()
    difference = difference.sort_values(by=[('seqNum', '')]).reset_index().drop([('index','')], axis=1)
    difference["xlabel"] = (difference[('seqNum', '')].astype(str) + "\n" + difference[('sequence','')].astype(str))
    difference = difference.drop([('seqNum', ''), ('sequence','')], axis=1)

    fig, ax = plt.subplots(figsize=(45, 6))

    ax.set_xticks(unidmeans.index[::10])
    ax.set_xticklabels(unidmeans["xlabel"][::10])

    for i in range(len(unidmeans)):
        if unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'YES':
            ax.axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='yellow',alpha=0.5)
        elif unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'YES':
            ax.axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='orange',alpha=0.3)

        bar_color1 = 'skyblue' if unidmeans.iloc[i, 0] >= 0 else 'lightgrey'
        bar_color2 = 'royalblue' if unidmeans.iloc[i, 1] >= 0 else 'lightgrey'

        ax.bar(unidmeans.index[i] - 0.225, unidmeans.iloc[i, 0], yerr=unidstdev.iloc[i, 0], capsize=2, width=0.45, color=bar_color1, align='center', label=f"{cond1_name}")
        ax.bar(unidmeans.index[i] + 0.225, unidmeans.iloc[i, 1], yerr=unidstdev.iloc[i, 1], capsize=2, width=0.45, color=bar_color2, align='center', label=f"{cond2_name}")

        ax.axhline(y=0.4, color="orange", linestyle="-", label="Medium Reactivity")
        ax.axhline(y=0.7, color="red", linestyle="-", label="High reactivity")
        ax.axhline(y=0.0, color="silver", linestyle="-")

    legend_elements = [
        Patch(facecolor="skyblue", edgecolor='grey',label=f"{cond1_name}"),
        Patch(facecolor="royalblue", edgecolor='grey', label=f"{cond2_name}"),
        Patch(facecolor="lightgrey", edgecolor='grey',label="Undetermined"),
        Patch(facecolor="yellow", edgecolor='grey',label=f"{cond2_name} sign. higher"),
        Patch(facecolor="orange", edgecolor='grey',label=f"{cond2_name} sign. lower"),
        Line2D([0], [0], color="red", label="High reactivity threshold"),
        Line2D([0], [0], color="orange", label="Medium reactivity threshold"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", ncol=14, fontsize=15)
    # ax.set_title(f"Nucleotides {unidmeans['xlabel'][0].split()[0]} - {unidmeans['xlabel'][-1].split()[0]}",loc='left', fontsize=20,loc='left', fontsize=20)
    ax.set_xlabel("Sequence", fontsize=20)
    ax.set_ylabel("Reactivity", fontsize=20)

    if unidstdev.iloc[:, 0].isna().all():
        unidstdev.iloc[:, 0] = 0
    if unidstdev.iloc[:, 1].isna().all():
        unidstdev.iloc[:, 1] = 0

    ax.set_xlim([unidmeans.index[0] - 1, unidmeans.index[-1] + 1])
    ax.set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 1]-unidstdev.iloc[:, 1]))*1.1,-0.15), \
        max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 1]+unidstdev.iloc[:, 1]))*1.2])
    ax2 = ax.twinx()
    ax2.set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 1]-unidstdev.iloc[:, 1]))*1.1,-0.15), \
        max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 1]+unidstdev.iloc[:, 1]))*1.2])
    ax.tick_params(labelsize=20)
    ax2.tick_params(labelsize=20)
    ax.set_title(title, fontsize=30)
    ax.grid(axis='y')
    # plt.suptitle(title, fontsize=30)
    plt.tight_layout()
    plt.savefig(output, format=format)

    fig, ax = plt.subplots(figsize=(45, 6))

    ax.set_xticks(unidmeans.index[::10])
    ax.set_xticklabels(unidmeans["xlabel"][::10])

    for i in range(len(unidmeans)):
        if unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'YES':
            ax.bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'red',edgecolor='red', align='center', label="Significant higher")
        elif unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'YES':
            ax.bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'blue',edgecolor='blue', align='center', label="Significant lower")
        elif unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'NO' and unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'NO':
            ax.bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'palegreen', edgecolor='forestgreen',align='center', label="Difference")
        else:
            ax.bar(unidmeans.index[i], -0.1, width=0.8, color = 'lightgrey', align='center', label="Undetermined")


    ax.axhline(y=0.0, color="silver", linestyle="-")

    legend_elements = [
        Patch(facecolor="palegreen", label="Difference"),
        Patch(facecolor="red", label="Significant higher"),
        Patch(facecolor="blue", label="Significant lower"),
        Patch(facecolor="lightgrey", label="Undetermined"),
    ]

    if difference.iloc[:, 0].isna().all():
        difference.iloc[:, 0] = 0

    ax.legend(handles=legend_elements, loc="upper left", ncol=14, fontsize=30)
    # ax.set_title(f"Nucleotides {unidmeans['xlabel'][0].split()[0]} - {unidmeans['xlabel'][-1].split()[0]}",loc='left', fontsize=20)
    ax.set_xlim([difference.index[0] - 1, difference.index[-1] + 1])
    ax.set_xlabel("Sequence", fontsize=20)
    ax.set_ylabel("Reactivity", fontsize=20)
    ax.set_ylim([min(np.nanmin(difference.iloc[:, 0])*1.1, -0.15), np.nanmax(difference.iloc[:, 0])*1.1])
    ax2 = ax.twinx()
    ax2.set_ylim([min(np.nanmin(difference.iloc[:, 0])*1.1, -0.15), np.nanmax(difference.iloc[:, 0])*1.1])
    ax.tick_params(labelsize=20)
    ax2.tick_params(labelsize=20)
    ax.set_title(dif_title, fontsize=30)
    ax.grid(axis='y')
    # plt.suptitle(dif_title, fontsize=30)
    plt.tight_layout()
    plt.savefig(diff_output, format=format)


def footprint_2D_plot(infile, higher, lower, outfile):
    vecs=[higher,lower]
    for i,vec in enumerate(vecs):
        pos=[i for i,x in vec.items() if x=='YES']
        vecs[i] = str(pos).replace(' ','')[1:-1]

    import os 
    conda_prefix = os.environ.get("CONDA_PREFIX")

    varna_cmd = ['java', '-jar', f'{conda_prefix}/lib/varna/VARNA.jar',
                 #'-sequenceDBN', sequence,
                 #'-structureDBN', structure,
                 '-i', infile,
                 '-o', outfile,
                 '-basesStyle1', 'fill=#ff0000',
                 '-basesStyle2', 'fill=#0000ff',]

    for i,vec in enumerate(vecs):
        if vec:
            varna_cmd.extend([f'-applyBasesStyle{i+1}on', vec])

    import subprocess
    subprocess.run(varna_cmd)

def footprint_main(
    path,
    # prefix,
    cond1_name: str = None,
    cond2_name: str = None,
    # deviation_type: str = "stdev",
    ttest_pvalue_thres=0.05,
    diff_thres=0.2,
    ratio_thres=0.2,
    output="fig.svg",
    plot: str = None,
    diff_plot="diff.svg",
    diff_plot_title="Footprint",
    plot_title="Footprint",
    plot_format="svg",
    # structure = None, ## dbn structure file
    # structure_plot=None,
):
    
    conditions = [(cond1_name, cond2_name), (cond2_name, cond1_name)]

    for cp1, cp2 in tqdm(conditions, total=len(conditions), desc="Writing files"):

        os.makedirs(f'{path}/{cp1}/comp_{cp1}_{cp2}', exist_ok=True)
        fnames = os.listdir(f'{path}/{cp1}')
        rnanames = [f.split('_')[-2] for f in fnames if f.split('_')[-1] == 'aggregated.tsv']
        for prefix in rnanames:
            os.makedirs(f'{path}/{cp1}/comp_{cp1}_{cp2}/{prefix}', exist_ok=True)
            footprint, footprint_csv = footprint_ttest(
                f'{path}/{cp1}/{cp1}_{prefix}_aggregated.tsv',
                cp1,
                f'{path}/{cp2}/{cp2}_{prefix}_aggregated.tsv',
                cp2,
                # deviation_type=deviation_type,
                ttest_pvalue_thres=0.05,
                diff_thres=0.2,
                ratio_thres=0.2,
            )
    
        # if output is not None:
            footprint_csv.to_csv(f'{path}/{cp1}/comp_{cp1}_{cp2}/{prefix}/{cp1}_{cp2}_{prefix}_footprint.tsv', sep="\t")
    
        # if plot is not None or diff_plot is not None:
            plot_reactivity(
                footprint,
                plot_title,
                diff_plot_title,
                plot_format,
                output=f'{path}/{cp1}/comp_{cp1}_{cp2}/{prefix}/{cp1}_{cp2}_{prefix}_footprint.svg',
                diff_output=f'{path}/{cp1}/comp_{cp1}_{cp2}/{prefix}/{cp1}_{cp2}_{prefix}_footprint_diff.svg',
                # deviation_type=deviation_type,
                cond1_name=cp1,
                cond2_name=cp2,
            )
    
        # if structure_plot is not None:
            # Write annotated structure plot {structure_plot} based on file {structure}
            structure = f'{path}/{cp1}/{cp1}_{prefix}.dbn'
    
            assert(structure is not None)
            footprint1 = footprint_csv.reset_index().set_index('seqNum')
            higher1 = footprint1['analysis']['significant_higher']
            lower1 = footprint1['analysis']['significant_lower']
            footprint_2D_plot(structure, higher1, lower1, outfile=f'{path}/{cp1}/comp_{cp1}_{cp2}/{prefix}/{cp1}_{cp2}_{prefix}_footprint_structure.svg')


def main():
    fire.Fire(footprint_main)


if __name__ == "__main__":
    main()



if __name__ == "__main__":
    main()
