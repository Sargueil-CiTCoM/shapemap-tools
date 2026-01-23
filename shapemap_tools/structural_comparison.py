import os
import re
import sys
import RNA
import pandas as pd
import numpy as np
import varnaapi
import fire
from collections import Counter



def structural_comparison(filename1, filename2):
    """
    Compare secondary structures.
    """

    def get_sequence_structure(filename):
        with open(f'{filename}', 'r') as file:
                content = file.read()
        lines = content.split('\n')
        sequence = lines[1]
        structure = lines[2]
        mfe_base_pair = []
        mfe_bp = list(RNA.ptable(structure))[1:]
        for j in range(len(mfe_bp)):
            if mfe_bp[j] != 0:
                mfe_base_pair.append([j,mfe_bp[j]-1])
        model = {'sequence': sequence, 'secondary_structure':structure, 'mfe_base_pair':mfe_base_pair}
        return model

    def model_vs_model(model1, model2):
        identical_bases = []
        different_bps = []
        for n, i in enumerate(model1['secondary_structure']):
            if i == model2['secondary_structure'][n] and model2['secondary_structure'][n] not in ['(',')']:
                identical_bases.append(1)
                different_bps.append(1)
            elif i in ['(',')'] and model2['secondary_structure'][n] in ['(',')']:
                bp = [bp for bp in model1['mfe_base_pair'] if bp[0] == n]
                if bp[0] in model2['mfe_base_pair']:
                    identical_bases.append(1)
                    different_bps.append(1)
                else:
                    identical_bases.append(0)
                    different_bps.append(-1)
            elif i == '.' and model2['secondary_structure'][n] == '-':
                identical_bases.append(1)
                different_bps.append(1)
            else:
                identical_bases.append(0)
                different_bps.append(0)
        counter_1 = dict(Counter(identical_bases))
        counter_2 = dict(Counter(different_bps))
        ratio_counter_1 = {('Different' if s == 0 else 'Identical' if s == 1 else 'Paired with a different partner'): counter_1[s]/len(model1['sequence']) for s in counter_1}
        ratio_counter_2 = {('Different' if s == 0 else 'Identical' if s == 1 else 'Paired with a different partner'): counter_2[s]/len(model1['sequence']) for s in counter_2} 
        model = {'sequence' : model1['sequence'], 'secondary_structure' : model1['secondary_structure'], 'identical_bases':identical_bases, 'different_bps':different_bps}
        proportion = {'proportion_identical_bases':ratio_counter_1, 'proportion_different_bps':ratio_counter_2}
        return model, proportion

    def varnaplot(model_info, color, path):
        conda_prefix = os.environ.get("CONDA_PREFIX")
        varnaapi.set_VARNA(f'{conda_prefix}/lib/varna/VARNA.jar')
        v = varnaapi.Structure(sequence=model_info['sequence'], structure=model_info['secondary_structure'])
        v.set_algorithm('radiate')
        v.update(bpStyle="lw", spaceBetweenBases=0.75, bpIncrement=1.3)
        v.add_colormap(values = color, vMin=-1, vMax=1, caption='', style={-1: '#fc0303',-0.9: '#ffffff', 0.9:'#ffffff',1:'#1d05f5'})
        tmpname = path
        v.savefig(tmpname)
            
    model1 = get_sequence_structure(filename1)
    model2 = get_sequence_structure(filename2)
    if model1['sequence'].upper().replace('T','U') != model2['sequence'].upper().replace('T','U') :
        print('\n*** Please verify the input DBN files, the sequences differ between the two DBN files. ***')

    elif len(model1['secondary_structure']) != len(model2['secondary_structure']):
        print('\n*** Please verify the input DBN files, the secondary structures in the two DBN files have different lengths. ***')
    
    else:
        output_name = input("Please enter the output folder name : ")
        while output_name == '' or output_name == ' ':
            output_name = input("Please enter the output folder name : ")
        if output_name != '' and output_name != ' ':
            output_name = output_name
            model_vs, proportion_vs = model_vs_model(model1, model2)

            os.makedirs(f'./comparaison_models/{output_name}', exist_ok=True)
            file_paths = {}
            for j in model_vs:
                if j not in ['sequence', 'secondary_structure']:
                    file_paths[f'{j}_file_path'] = f'./comparaison_models/{output_name}/{j}.txt'
                    varnaplot_path = f'./comparaison_models/{output_name}/{j}.varna'
                    varnaplot(model_vs, model_vs[j], varnaplot_path)
            for pp in proportion_vs:
                with open(f'./comparaison_models/{output_name}/{pp}.txt', "w", encoding="utf-8") as f:
                    f.write('Proportion :\n')
                    for key, value in proportion_vs[pp].items():
                        f.write(f"{key}: {value}\n")


def main():
    fire.Fire(structural_comparison)


if __name__ == "__main__":
    main()
