#!/bin/bash


# The LIBS variable contains the following columns : 
# exp_name,Modified_id,Unmodified_id,fasta-sequence,


LIBS="\
    SAMAP_rep1-nosam_Mg5,S3C,S1A,S1E,../data/sequences/2022-05-13-aptmer-full_riboswitch-sequence-for-shapemapper-no-t7.fasta   \
    SAMAP_rep1-sam01_Mg5,S3G,S1A,S1E,../data/sequences/2022-05-13-aptmer-full_riboswitch-sequence-for-shapemapper-no-t7.fasta   \
    SAMAP_rep1-sam1_Mg5,S3I,S1A,S1E,../data/sequences/2022-05-13-aptmer-full_riboswitch-sequence-for-shapemapper-no-t7.fasta    \
    "
# Where to find the input data
DATA_FOLDER1="../data/2022-05-13-decrypted-aptamers-riboswitchs"
DATA_FOLDER2="../data/2022-06-29-decrypted_aptamers_riboswitchs_retry"
DATA_FOLDER3="../data/2022-08-12-VIH-Ribo-aptamers-PDB"

# Number of data folder
NDATA_FOLDER=3



ARN_FILE=tags.tsv
SPLIT_ARN=1


# Where to store results
OUTPUT_PATH="../data/results/"

# Arguments to shape mapper
SHAPEMAPPER_ARGS="--nproc 48 --indiv-norm --min-depth 5000"

# Name of the output
SHAPEMAPPEROUT_NAME="2022-08-12-aptamer-ribo_Rep3_5000depth"

# Wait for a program before starting treatement
WAIT_PID=

source run.inc.sh
