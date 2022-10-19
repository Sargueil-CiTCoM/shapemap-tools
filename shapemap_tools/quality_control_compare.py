#! /usr/bin/env python3

import fire
import subprocess
import os

base_path = os.path.dirname(__file__)


def report(data_path1, data_path2, mode="condition", dev=False):
    path = os.path.join(base_path, "qualitycontrol_compare_report.ipynb")
    env = os.environ.copy()
    env["DATA_PATH1"] = os.path.abspath(data_path1)
    env["DATA_PATH2"] = os.path.abspath(data_path2)
    print(data_path1)
    print(data_path2)
    env["MODE"] = mode
    if dev:
        subprocess.Popen(["jupyter-notebook", path], env=env)
    else:
        subprocess.Popen(["voila", path], env=env)


def main():
    """TODO: Docstring for main.
    Returns
    -------
    TODO

    """
    fire.Fire(report)


if __name__ == "__main__":
    main()
