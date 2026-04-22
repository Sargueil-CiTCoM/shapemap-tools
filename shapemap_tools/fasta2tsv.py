import pandas as pd
import fire
import os
import re
import sys
import os


def fasta2tsv(input_fasta, output_tsv):
    """
    Convert a FASTA file to a TSV file.

    Args:
        input_fasta: path to input FASTA file
        output_tsv: path to output TSV file
    """

    sequences = {}
    rnaname = None
    with open(input_fasta, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                rnaname = line[1:]
                sequences[rnaname] = ''
            elif rnaname:
                sequences[rnaname] += line
    
    df = pd.DataFrame(sequences.items(), columns=['name', 'sequence'])
    df.to_csv(output_tsv, index=False, sep='\t')



def main():
    fire.Fire(fasta2tsv)


if __name__ == "__main__":
    main()


