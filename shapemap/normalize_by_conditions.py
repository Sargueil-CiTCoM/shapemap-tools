#!/usr/bin/env python3
# coding: utf-8

# In[1]:


import os
import fire
import subprocess as sp
import glob

# import multiprocessing as mp


shapemapper_path = "../shapemapper2/internals/bin/"
norm_script_path = os.path.join(shapemapper_path, "normalize_profiles.py")
tab_to_shape_path = os.path.join(shapemapper_path, "tab_to_shape.py")
render_figures_path = os.path.join(shapemapper_path, "render_figures.py")

base_path = "{path}/{condition}/"
profile_pattern = base_path + "{condition}_{seqid}_profile.txt"
shape_pattern = base_path + "{condition}_{seqid}.shape"
map_pattern = base_path + "{condition}_{seqid}.map"
histo_pattern = base_path + "{condition}_{seqid}_histograms.pdf"
plot_pattern = base_path + "{condition}_{seqid}_profiles.pdf"
title_pattern = "{seqid} {condition}"


def runNormalizeShapemapper(tonorm, normout):
    cmd = ["python", norm_script_path, "--tonorm"] + tonorm + ["--normout"] + normout
    # [print(arg,end=" ") for arg in cmd]
    try:
        sp.run(cmd)
    except Exception:
        pass


def runTabToShape(profile_file, map_file, shape_file):
    cmd = [
        "python",
        tab_to_shape_path,
        "--infile",
        profile_file,
        "--map",
        map_file,
        "--shape",
        shape_file,
    ]
    try:
        sp.run(cmd)
    except Exception:
        pass


def runRenderFigures(
    profile_file, min_depth=5000, title="", plot_file=None, histo_file=None
):
    cmd = [
        "python",
        render_figures_path,
        "--infile",
        profile_file,
        "--mindepth",
        str(min_depth),
        "--title",
        title,
        "--hist",
        histo_file,
        "--plot",
        plot_file,
    ]
    try:
        sp.run(cmd)
    except Exception:
        pass


def normalize_through_all(path, outputpath, sequences, conditions):
    tonorm = [
        profile_pattern.format(path=path, condition=condition, seqid=seqid)
        for condition in conditions
        for seqid in sequences
        if os.path.exists(
            profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    normout = [
        profile_pattern.format(path=outputpath, condition=condition, seqid=seqid)
        for condition in conditions
        for seqid in sequences
        if os.path.exists(
            profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    runNormalizeShapemapper(tonorm, normout)

    for cond in conditions:
        for seqid in sequences:
            print(f"Start {cond} ({seqid})")
            runTabToShape(
                profile_pattern.format(path=outputpath, condition=cond, seqid=seqid),
                map_file=map_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
                shape_file=shape_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
            )
            runRenderFigures(
                profile_pattern.format(
                    path=outputpath,
                    condition=cond,
                    seqid=seqid,
                ),
                title=title_pattern.format(condition=cond, seqid=seqid),
                plot_file=plot_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
                histo_file=histo_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
            )


def normalize_through_conditions(path, outputpath, seqid, conditions):

    tonorm = [
        profile_pattern.format(path=path, condition=condition, seqid=seqid)
        for condition in conditions
        if os.path.exists(
            profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    normout = [
        profile_pattern.format(path=outputpath, condition=condition, seqid=seqid)
        for condition in conditions
        if os.path.exists(
            profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    print(f"Start {seqid}")
    runNormalizeShapemapper(tonorm, normout)
    for cond in conditions:
        print(f"Start {cond} ({seqid})")
        runTabToShape(
            profile_pattern.format(path=outputpath, condition=cond, seqid=seqid),
            map_file=map_pattern.format(path=outputpath, condition=cond, seqid=seqid),
            shape_file=shape_pattern.format(
                path=outputpath, condition=cond, seqid=seqid
            ),
        )
        runRenderFigures(
            profile_pattern.format(
                path=outputpath,
                condition=cond,
                seqid=seqid,
            ),
            title=title_pattern.format(condition=cond, seqid=seqid),
            plot_file=plot_pattern.format(path=outputpath, condition=cond, seqid=seqid),
            histo_file=histo_pattern.format(
                path=outputpath, condition=cond, seqid=seqid
            ),
        )


def main(globpath, outputpath, mode="conditions"):

    os.makedirs(outputpath, exist_ok=True)
    conditions = [os.path.basename(path) for path in glob.glob(f"{globpath}/*")]

    sequences = set()
    for cond in conditions:
        os.makedirs(os.path.join(outputpath, cond), exist_ok=True)
        sequences = sequences.union(
            {
                os.path.splitext(os.path.basename(path))[0][(len(cond) + 1) :]
                for path in glob.glob(
                    shape_pattern.format(path=globpath, condition=cond, seqid="*")
                )
            }
        )

    if mode == "conditions":
        for seqid in sequences:
            print(f"Prep {seqid}")
            args = (globpath, outputpath, seqid, conditions)
            normalize_through_conditions(*args)
            # pool.apply_async(normalize_through_conditions, args)
    elif mode == "all":
        args = (globpath, outputpath, sequences, conditions)
        normalize_through_all(*args)
        # pool.apply_async(normalize_through_conditions, args)
    else:
        raise Exception(f"Invalid mode {mode} Your must choose between: all, conditions")


#    #pool.close()
#    #pool.join()
#

if __name__ == "__main__":

    fire.Fire(main)
