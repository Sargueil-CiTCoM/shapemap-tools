import subprocess as sp
import fire
import utils

def runIPANEMAP(configfile):
    cmd = ["ipanemap", "--config", configfile]

    # [print(arg,end=" ") for arg in cmd]
    sp.run(cmd)


def gen_structures():
    """TODO: Docstring for gen_structures.
    Returns
    -------
    TODO

    """
    pass


def main():
    fire.Fire(gen_structures)


if __name__ == "__main__":
    main()
