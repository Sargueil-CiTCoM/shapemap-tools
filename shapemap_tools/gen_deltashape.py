#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import re
import os
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import SVG, display
from scipy import sparse
import skbio as sb
import subprocess as sp
import sys
sys.path.append("..")



def runDeltaShape(profile1, profile2, output, title, smoothing=1, find_site=(2,3), low_thres=-999, invalid_low_thres=-999):
    cmd = ["python", deltashape_path, profile1, profile2, "-o", output, "--title", title, "-f", f"{find_site[0]},{find_site[1]}", "--pad", str(smoothing),"--low-thres", str(low_thres), "--invalid-low-thres", str(invalid_low_thres), "--noshow", "--all", "--svg", "--png", ]
    #[print(arg,end=" ") for arg in cmd]
    sp.run(cmd)

def genDelta(config, dsname1, dsname2, smoothing=2, find_site=(2,3), low_thres=-999, invalid_low_thres=-999):

    path1= config.shapem_cond_path(dsname1)
    path2 = config.shapem_cond_path(dsname2)
    ids = config.aptamers()
    
    dspath = config.deltashape_cond_path(dsname1, dsname2)
    
    os.makedirs(f'{dspath}', exist_ok=True)
    os.makedirs(f'{dspath}/svg', exist_ok=True)
    os.makedirs(f'{dspath}/png', exist_ok=True)
    for seqid in ids:
        p1 = f"{path1}/{dsname1}_{seqid}.map"
        p2 = f"{path2}/{dsname2}_{seqid}.map"
        output = f"{dspath}/deltashape_{seqid}_{dsname1}_{dsname2}"
        title = f"{seqid} - Delta between {dsname1} and {dsname2}"
        #try:
        if True:
            runDeltaShape(p1, p2, f"{output}.tsv", title, smoothing=smoothing, find_site=find_site, low_thres=low_thres, invalid_low_thres=invalid_low_thres)
            sp.run(["mv", f"{output}.svg", f"{dspath}/svg"])
            sp.run(["mv", f"{output}.png", f"{dspath}/png"])
            print(title)
            #show_svg(f"{output}.svg")
        
        #except:
            #print("Unexpected error:", sys.exc_info()[0])
        #    print(f"Invalid data for {seqid}")
