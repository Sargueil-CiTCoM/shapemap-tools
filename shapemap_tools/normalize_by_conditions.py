#!/usr/bin/env python3
# coding: utf-8

# In[1]:


import os
import fire
import subprocess as sp
import glob
from . import utils

# import multiprocessing as mp


shapemapper_bin_path = "{shapemapper_path}/internals/bin/"
norm_script_path = "{shapemapper_bin_path}/normalize_profiles.py"
tab_to_shape_path = "{shapemapper_bin_path}/tab_to_shape.py"
render_figures_path = "{shapemapper_bin_path}/render_figures.py"


def runNormalizeShapemapper(tonorm, normout):
    cmd = (
        [
            "python",
            norm_script_path.format(shapemapper_bin_path=shapemapper_bin_path),
            "--tonorm",
        ]
        + tonorm
        + ["--normout"]
        + normout
    )
    # [print(arg,end=" ") for arg in cmd]
    try:
        return sp.run(cmd, capture_output=True, text=True).returncode
    except Exception as e:
        raise e
#        pass


def runTabToShape(profile_file, map_file, shape_file):
    cmd = [
        "python",
        tab_to_shape_path.format(shapemapper_bin_path=shapemapper_bin_path),
        "--infile",
        profile_file,
        "--map",
        map_file,
        "--shape",
        shape_file,
    ]
    try:
        return sp.run(cmd, capture_output=True, text=True).returncode
    except Exception as e:
        raise e
#        pass


def runRenderFigures(
    profile_file, min_depth=5000, title="", plot_file=None, histo_file=None
):
    cmd = [
        "python",
        render_figures_path.format(shapemapper_bin_path=shapemapper_bin_path),
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
        return sp.run(cmd, capture_output=True, text=True).returncode
    except Exception as e:
        raise e
#        pass


def normalize_through_all(path, outputpath, sequences, conditions):
    tonorm = [
        utils.profile_pattern.format(path=path, condition=condition, seqid=seqid)
        for condition in conditions
        for seqid in sequences
        if os.path.exists(
            utils.profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    normout = [
        utils.profile_pattern.format(path=outputpath, condition=condition, seqid=seqid)
        for condition in conditions
        for seqid in sequences
        if os.path.exists(
            utils.profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    runNormalizeShapemapper(tonorm, normout)

    for cond in conditions:
        for seqid in sequences:
            print(f"Start {cond} ({seqid})")
            runTabToShape(
                utils.profile_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
                map_file=utils.map_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
                shape_file=utils.shape_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
            )
            runRenderFigures(
                utils.profile_pattern.format(
                    path=outputpath,
                    condition=cond,
                    seqid=seqid,
                ),
                title=utils.title_pattern.format(condition=cond, seqid=seqid),
                plot_file=utils.plot_profiles_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
                histo_file=utils.plot_histo_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
            )


def normalize_through_conditions(path, outputpath, seqid, conditions):

    tonorm = [
        utils.profile_pattern.format(path=path, condition=condition, seqid=seqid)
        for condition in conditions
        if os.path.exists(
            utils.profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    normout = [
        utils.profile_pattern.format(path=outputpath, condition=condition, seqid=seqid)
        for condition in conditions
        if os.path.exists(
            utils.profile_pattern.format(path=path, condition=condition, seqid=seqid)
        )
    ]

    print(f"Start {seqid}")
    retcode = runNormalizeShapemapper(tonorm, normout)
    if retcode == 0:
        for cond in conditions:
            print(f"Start {cond} ({seqid})")
            tts_retcode = runTabToShape(
                utils.profile_pattern.format(path=outputpath, condition=cond, seqid=seqid),
                map_file=utils.map_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
                shape_file=utils.shape_pattern.format(
                    path=outputpath, condition=cond, seqid=seqid
                ),
            )
            if tts_retcode == 0:
                runRenderFigures(
                    utils.profile_pattern.format(
                        path=outputpath,
                        condition=cond,
                        seqid=seqid,
                    ),
                    title=utils.title_pattern.format(condition=cond, seqid=seqid),
                    plot_file=utils.plot_profiles_pattern.format(
                        path=outputpath, condition=cond, seqid=seqid
                    ),
                    histo_file=utils.plot_histo_pattern.format(
                        path=outputpath, condition=cond, seqid=seqid
                    ),
                )
            else:
                print(f"Shape conversion error for sequence {seqid} - {cond}")
    else:
        print(f"Normalization error for sequence {seqid}")


def main(globpath, outputpath, mode="conditions", shapemapper_path="shapemapper2"):
    """main

    normalize by condition
    Parameters
    ----------
    globpath :
        globpath
    outputpath :
        outputpath
    mode :
        mode
    shapemapper_path :
        shapemapper_path
    """

    global shapemapper_bin_path
    shapemapper_bin_path = shapemapper_bin_path.format(
        shapemapper_path=shapemapper_path
    )
    os.makedirs(outputpath, exist_ok=True)
    conditions = [os.path.basename(path) for path in glob.glob(f"{globpath}/*")]

    sequences = set()
    for cond in conditions:
        os.makedirs(os.path.join(outputpath, cond), exist_ok=True)
    sequences = utils.get_sequences(globpath, conditions)

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
        raise Exception(
            f"Invalid mode {mode} Your must choose between: all, conditions"
        )


#    #pool.close()
#    #pool.join()
#
def main_wrapper():
    fire.Fire(main)


if __name__ == "__main__":
    main_wrapper()
