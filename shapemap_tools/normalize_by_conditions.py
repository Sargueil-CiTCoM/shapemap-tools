# /usr/bin/env python3
# coding: utf-8

# In[1]:


import os
import fire
import subprocess as sp
import glob
from . import utils
from tqdm import tqdm

# import multiprocessing as mp
import logging

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
        res = sp.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            logging.error(" ".join(cmd))
            logging.error(res.stderr)
            logging.error(res.stdout)
        return res.returncode
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
        res = sp.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            logging.error(" ".join(cmd))
            logging.error(res.stderr)
            logging.error(res.stdout)
        return res.returncode
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
        res = sp.run(cmd, capture_output=True, text=True)
        if True:  # res.returncode != 0:
            logging.error(" ".join(cmd))
            logging.error(res.stderr)
            logging.error(res.stdout)

    except Exception as e:
        raise e

    try:
        res1 = sp.run(
            ["pdf2svg", plot_file, os.path.splitext(plot_file)[0] + ".svg"],
            capture_output=True,
            text=True,
        )
        res2 = sp.run(
            ["pdf2svg", histo_file, os.path.splitext(histo_file)[0] + ".svg"],
            capture_output=True,
            text=True,
        )

        if True:  # res.returncode != 0:
            logging.warning(res1.stderr)
            logging.warning(res1.stdout)
            logging.warning(res2.stderr)
            logging.warning(res2.stdout)

    except Exception as e:
        print(e)

    return res.returncode


def normalize_through_all(path, outputpath, sequences, conditions, normalize=True, render_figure=True):
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
    if normalize:
        runNormalizeShapemapper(tonorm, normout)

    for cond in tqdm(
        conditions,
        total=len(conditions),
        desc="conditions for {seqid}",
        position=1,
        leave=False,
    ):
        for seqid in sequences:
            print(f"Start {cond} ({seqid})")
            if normalize:
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
            if render_figure:
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


def normalize_through_conditions(path, outputpath, seqid, conditions, normalize=True, render_figure=True):

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
    if normalize:
        retcode = runNormalizeShapemapper(tonorm, normout)
    retcode = 0
    if retcode == 0:
        for cond in tqdm(
            conditions,
            total=len(conditions),
            desc=f"conditions for {seqid}",
            position=1,
            leave=False,
        ):

            if normalize:
                tts_retcode = runTabToShape(
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
            tts_retcode = 0
            if render_figure:
                if tts_retcode == 0:
                    rrf_retcode = runRenderFigures(
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
                    if rrf_retcode != 0:
                        print(f"Rendering figure failed {seqid} - {cond}")

                else:
                    print(f"Shape conversion error for sequence {seqid} - {cond}")
    else:
        print(f"Normalization error for sequence {seqid}")


def split_conds_by_rep(conditions):
    cond_by_reps = {}

    for cond in conditions:
        rep = cond.split("_")[-1]
        if rep not in cond_by_reps:
            cond_by_reps[rep] = []
        cond_by_reps[rep].append(cond)

    return cond_by_reps


def main(
    globpath,
    outputpath,
    mode="conditions",
    shapemapper_path="shapemapper2",
    condition_prefix="",
    log="normalize_by_cond.log",
    normalize=True,
    render_figure=True
):
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

    logging.basicConfig(level=logging.DEBUG, filename=log)
    global shapemapper_bin_path
    shapemapper_bin_path = shapemapper_bin_path.format(
        shapemapper_path=shapemapper_path
    )
    os.makedirs(outputpath, exist_ok=True)
    conditions = [
        os.path.basename(path) for path in glob.glob(f"{globpath}/{condition_prefix}*")
    ]

    sequences = set()
    for cond in tqdm(
        conditions,
        total=len(conditions),
        desc="Preparing paths",
        position=0,
        leave=False,
    ):
        os.makedirs(os.path.join(outputpath, cond), exist_ok=True)
    sequences = utils.get_sequences(globpath, conditions)

    if mode == "conditions":
        for seqid in tqdm(
            sequences, total=len(sequences), desc="Normalizing", position=0
        ):
            args = (globpath, outputpath, seqid, conditions, normalize, render_figure)
            normalize_through_conditions(*args)
            # pool.apply_async(normalize_through_conditions, args)
    elif mode == "conditions-reps":
        for seqid in tqdm(
            sequences, total=len(sequences), desc="Normalizing", position=0
        ):

            cond_by_reps = split_conds_by_rep(conditions)
            for name, cbr in cond_by_reps.items():
                args = (globpath, outputpath, seqid, cbr, normalize, render_figure)
                normalize_through_conditions(*args)
    elif mode == "all":
        args = (globpath, outputpath, sequences, conditions, normalize, render_figure)
        normalize_through_all(*args)
        # pool.apply_async(normalize_through_conditions, args)
    else:
        raise Exception(
            f"Invalid mode {mode} Your must choose between: all, conditions, conditions-reps"
        )


#    #pool.close()
#    #pool.join()
#
def main_wrapper():
    fire.Fire(main)


if __name__ == "__main__":
    main_wrapper()
