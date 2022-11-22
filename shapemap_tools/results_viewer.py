#! /usr/bin/env python3

import fire
import subprocess
import os

base_path = os.path.dirname(__file__)


def report(data_path, conditions_prefix, dev=False):
    path = os.path.join(base_path, "results_viewer.ipynb")
    env = os.environ.copy()
    env["DATA_PATH"] = os.path.abspath(data_path)
    env["CONDITIONS_PREFIX"] = conditions_prefix
    if dev:
        subprocess.run(["jupyter-notebook", path], env=env)
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
