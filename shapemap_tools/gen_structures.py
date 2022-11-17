import subprocess as sp
import multiprocessing as mp
import os
import fire
import glob
import parse
from . import fasta

varna_path = os.path.join (os.path.dirname(__file__),"VARNAcmd.jar")


def run_varna_thread(config, condition, rna, sequence, shape, outputdir):
    dbn = os.path.join(outputdir, condition, f"{condition}_{rna}.dbn")
    if os.path.exists(dbn):
        runVARNA(
            dbn,
            shape,
            f"{condition} - {rna}",
            os.path.join(outputdir, condition, f"{condition}_{rna}.varna"),
        )
        runVARNA(
            dbn,
            shape,
            f"{condition} - {rna}",
            os.path.join(outputdir, condition, f"{condition}_{rna}.svg"),
        )
    else:
        print(f"{dbn} not created - cannot launch varna")


def runVARNA(struct_file, shape_file, title, output_file, resolution=1.0):
    colormap = (
        "$-10.00:#CCCCCC,$-0.3001:#999999,$-0.30:#FFFFFF,$0.39999:#FFFFFF"
        ",$0.40:#FFFF47,$0.8499:#FFFF47,$0.85:#FF0000,10:#FF0000"
    )

    cmd = [
        "java",
        "-jar",
        varna_path,
        "-i",
        struct_file,
        "-o",
        output_file,
        "-colorMapStyle",
        colormap,
        "-colorMap",
        shape_file,
        "-title",
        title,
        "-resolution",
        str(resolution),
    ]
    sp.run(cmd)


def run_thread(config, condition, rna, sequence, shape, outputdir):
    dbn = os.path.join(outputdir, condition, f"{condition}_{rna}.dbn")
    run_ipanemap(
        config,
        condition,
        sequence,
        shape,
        dbn,
        os.path.join(outputdir, condition, "tmp"),
    )


def run_ipanemap(config, condition, sequence, shape, dbn, tmp_outputdir):
    cmd = [
        "ipanemap",
        "--config",
        config,
        "-s",
        sequence,
        "--dbn-pattern",
        dbn,
        "--def-conditions",
        f"{condition}::{shape}:",
        "--conditions",
        condition,
        "--out-dir",
        tmp_outputdir,
    ]
    print(" ".join(cmd))
    print(f"Starting {condition}")
    sp.run(cmd)
    print(f"End {condition}")


def gen_structures(
    input_dir,
    output_dir,
    sequence,
    configfile=None,
):
    if configfile is None:
        configfile = os.path.join(
            os.path.dirname(__file__), "template", "ipanemap-config.yaml"
        )

    sequences = dict(
        [(name, seq.replace("T", "U")) for name, seq in fasta.fasta_iter(sequence)]
    )
    os.makedirs(f"{output_dir}/sequences", exist_ok=True)
    for name, seq in sequences.items():
        with open(f"{output_dir}/sequences/{name}.fasta", "w") as out:
            out.write(f"> {name}\n")
            out.write(seq)

    shape_files_path = glob.glob(f"{input_dir}/**/*.shape")
    conditions = set()

    for file in shape_files_path:
        conditions.add(os.path.basename(os.path.dirname(file)))

    for cond in conditions:
        os.makedirs(f"{output_dir}/{cond}/tmp", exist_ok=True)
    # pool = mp.Pool(24)
    for file in shape_files_path:
        base_name = os.path.basename(parse.parse("{name}.shape", file).named["name"])
        condname = os.path.basename(os.path.dirname(file))
        rna_name = base_name.split("_")[-1]
        kwargs = {
            "config": configfile,
            "condition": condname,
            "rna": rna_name,
            "sequence": f"{output_dir}/sequences/{rna_name}.fasta",
            "shape": file,
            "outputdir": f"{output_dir}",
        }
        # run_thread(
        #    **kwargs,
        # )
        run_thread(**kwargs)
        run_varna_thread(**kwargs)
        # pool.apply_async(
        #    run_thread,
        #    kwds=kwargs,
        # )
    # pool.close()
    # pool.join()


def main():
    fire.Fire(gen_structures)


if __name__ == "__main__":
    main()
