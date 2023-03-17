#!/bin/bash
SHAPEMAPPER_BIN="../../repos/shapemapper2"
DECRYPTED_FOLDER=../../repos/decrypted/data
DEMUX_FOLDER=demux_fastq
RAW_FOLDER=../raw_fastq
SAMAP_PLUS_PDB_INFO=${DECRYPTED_FOLDER}/1_tagged_sequences/2021-12-16-tagged-aptamer-sequences-natural_and_synthetic_with_source_family-corrected.tsv
RIBOSAM_INFO=${DECRYPTED_FOLDER}/1_tagged_sequences/2021-09-15-tagged-riboswitch-sam-natural.tsv
RIBOSAM_PDB_INFO=${DECRYPTED_FOLDER}/1_tagged_sequences/2022-01-12-tagged-RIBOSAMN-PDB.tsv
SAMAP_PDB_INFO=${DECRYPTED_FOLDER}/1_tagged_sequences/2022-05-13-tagged-SAMAP-PDB-modified.tsv

SAMAP_PLUS_PDB_INFO=${DECRYPTED_FOLDER}/1_tagged_sequences/2022-10-27-tagged-aptamer-sequences-natural_artifical_pdb.tsv
RIBOSAM_PLUS_PDB_INFO=${DECRYPTED_FOLDER}/1_tagged_sequences/2022-10-27-tagged-riboswitch-sam-natural_pdb.tsv

SAMAP_PLUS_PDB_FASTA=results/SAMAP/sequences/2022-10-27-tagged-aptamer-sequences-natural_artifical_pdb.fasta
RIBOSAM_PLUS_PDB_FASTA=results/RIBOSAM/sequences/2022-10-27-tagged-riboswitch-sam-natural_pdb.fasta


WAIT_FOR_SCRIPT=NONE #"parnassus.sh"

function is_up()
{
  NPROC=`ps aux | grep "$1" | grep -v grep | wc -l`
  test ${NPROC} -ne 0
}


while is_up $WAIT_FOR_SCRIPT; do
	echo \"$WAIT_FOR_SCRIPT\" is still running
	echo sleeping for 60s
	sleep 60;
done;

## DEMULTIPLEXING
# SAMAP
set -e
set -x


# rep0
#

echo "-- DEMUX -- SAMAP --"
shpm_demultiplex $SAMAP_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2021-11-23-decrypted-SAMAP-rep0 \
    --folder ${RAW_FOLDER}/2021-11-23-riboswitch-ngs-results;

    # rep1 - rep2

shpm_demultiplex $SAMAP_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-05-13-decrypted-SAMAP-rep1-1 \
    --folder ${RAW_FOLDER}/2022-05-13-decrypted-aptamers-riboswitchs

shpm_demultiplex $SAMAP_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-06-29-decrypted-SAMAP-rep1-2 \
    --folder ${RAW_FOLDER}/2022-06-29-decrypted_aptamers_riboswitchs_retry

# rep3
shpm_demultiplex $SAMAP_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-08-10-decrypted-SAMAP-rep3 \
    --folder ${RAW_FOLDER}/2022-08-10-aptamer

echo "-- DEMUX -- RIBOSAM --"
# RIBOSAM

#rep 0
shpm_demultiplex $RIBOSAM_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-05-13-decrypted-RIBOSAM-rep1-1 \
    --folder ${RAW_FOLDER}/2022-05-13-decrypted-aptamers-riboswitchs

shpm_demultiplex $RIBOSAM_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-06-29-decrypted-RIBOSAM-rep1-2 \
    --folder ${RAW_FOLDER}/2022-06-29-decrypted_aptamers_riboswitchs_retry

# rep1
shpm_demultiplex $RIBOSAM_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-08-10_12-decrypted-RIBOSAM-rep0 \
    --folder ${RAW_FOLDER}/2022-08-10_12_ribosam-rep0
        ## SHAPEMAPPER

 rep3
shpm_demultiplex $RIBOSAM_PLUS_PDB_INFO --cores 48 --max-errors 1 \
    --output_path demux_fastq/2022-08-10-decrypted-RIBOSAM-rep3 \
    --folder ${RAW_FOLDER}/2022-08-10-aptamer

echo "-- SHAPEMAPPER -- SAMAP --"
shpm_launch_shapemapper \
    config/decrypted-SAMAP_config.yaml \
    config/decrypted-SAMAP_samples.tsv  --shapemapper-path $SHAPEMAPPER_BIN/shapemapper --dnerase=True --verbose=True --nthreads=48


echo "-- SHAPEMAPPER -- RIBOSAM --"
shpm_launch_shapemapper \
        config/decrypted-RIBOSAM_config.yaml \
        config/decrypted-RIBOSAM_samples.tsv --shapemapper-path $SHAPEMAPPER_BIN/shapemapper --dnerase=True --nthreads=48


echo "-- NORM_ALL_REPS -- SAMAP --"
shpm_norm_by_conditions results/SAMAP results/SAMAP-norm-conds --shapemapper-path $SHAPEMAPPER_BIN --condition-prefix SAMAP --mode conditions --render-figure=False --nthreads=48

echo "-- NORM_ALL_REPS -- RIBOSAM --"
shpm_norm_by_conditions results/RIBOSAM results/RIBOSAM-norm-conds --shapemapper-path $SHAPEMAPPER_BIN --condition-prefix RIBOSAM --mode conditions --render-figure=False --nthreads=48

echo "-- NORM_REPS -- SAMAP --"
shpm_norm_by_conditions results/SAMAP results/SAMAP-norm-conds-by-reps --shapemapper-path $SHAPEMAPPER_BIN --condition-prefix SAMAP --mode conditions-reps --render-figure=False --nthreads=48
echo "-- NORM_REPS -- RIBOSAM --"
shpm_norm_by_conditions results/RIBOSAM results/RIBOSAM-norm-conds-by-reps --shapemapper-path $SHAPEMAPPER_BIN --condition-prefix RIBOSAM --mode conditions-reps --render-figure=False --nthreads=48


echo "-- AGGREGATE -- SAMAP --"
shpm_agg_replicates --input_path results/SAMAP --output_path results/SAMAP-aggregated --config_path config/decrypted-SAMAP_config.yaml


echo "-- AGGREGATE -- RIBOSAM --"
shpm_agg_replicates --input_path results/RIBOSAM --output_path results/RIBOSAM-aggregated --config_path config/decrypted-RIBOSAM_config.yaml

echo "-- AGGREGATE -- SAMAP_NORMED_ALL_REPS --"
shpm_agg_replicates --input_path results/SAMAP-norm-conds --output_path results/SAMAP-norm-conds-aggregated --config_path config/decrypted-SAMAP_config.yaml

echo "-- AGGREGATE -- RIBOSAM_NORMED_ALL_REPS --"
shpm_agg_replicates --input_path results/RIBOSAM-norm-conds --output_path results/RIBOSAM-norm-conds-aggregated --config_path config/decrypted-RIBOSAM_config.yaml

echo "-- AGGREGATE -- SAMAP_NORMED_REPS --"
shpm_agg_replicates --input_path results/SAMAP-norm-conds-by-reps --output_path results/SAMAP-norm-conds-by-reps-aggregated --config_path config/decrypted-SAMAP_config.yaml
echo "-- AGGREGATE -- RIBOSAM_NORMED_REPS --"
shpm_agg_replicates --input_path results/RIBOSAM-norm-conds-by-reps --output_path results/RIBOSAM-norm-conds-by-reps-aggregated --config_path config/decrypted-RIBOSAM_config.yaml

echo "-- PLOT_AGGREGATE -- SAMAP --"
shpm_plot_aggregate results/SAMAP-aggregated --nthreads=48 #--dnerase
echo "-- PLOT_AGGREGATE -- RIBOSAM --"
shpm_plot_aggregate results/RIBOSAM-aggregated --dnerase

echo "-- PLOT_AGGREGATE -- SAMAP_NORMED_ALL_REPS --"
shpm_plot_aggregate results/SAMAP-norm-conds-aggregated --nthreads=48 #--dnerase
echo "-- PLOT_AGGREGATE -- RIBOSAM_NORMED_ALL_REPS --"
shpm_plot_aggregate results/RIBOSAM-norm-conds-aggregated --dnerase

echo "-- PLOT_AGGREGATE -- SAMAP_NORMED_REPS --"
shpm_plot_aggregate results/SAMAP-norm-conds-by-reps-aggregated --nthreads=48 #--dnerase
echo "-- PLOT_AGGREGATE -- RIBOSAM_NORMED_REPS --"
shpm_plot_aggregate results/RIBOSAM-norm-conds-by-reps--aggregated --dnerase


echo "-- GENSTRUCT -- SAMAP --"
shpm_gen_structures results/SAMAP-aggregated results/SAMAP-aggregated $SAMAP_PLUS_PDB_FASTA --nthreads=24
echo "-- GENSTRUCT -- RIBOSAM --"
shpm_gen_structures results/RIBOSAM-aggregated results/RIBOSAM-aggregated $RIBOSAM_PLUS_PDB_FASTA

echo "-- GEN_STRUCT -- SAMAP_NORMED_ALL_REPS --"
shpm_gen_structures results/SAMAP-norm-conds-aggregated results/SAMAP-norm-conds-aggregated $SAMAP_PLUS_PDB_FASTA --nthreads=24
echo "-- GEN_STRUCT -- RIBOSAM_NORMED_ALL_REPS --"
shpm_gen_structures results/RIBOSAM-norm-conds-aggregated results/RIBOSAM-norm-conds-aggregated $RIBOSAM_PLUS_PDB_FASTA

echo "-- GENSTRUCT -- SAMAP_NORMED_REPS --"
shpm_gen_structures results/SAMAP-norm-conds-by-reps-aggregated results/SAMAP-norm-conds-by-reps-aggregated $SAMAP_PLUS_PDB_FASTA --nthreads=24
echo "-- GENSTRUCT -- RIBOSAM_NORMED_REPS --"
shpm_gen_structures results/RIBOSAM-norm-conds-by-reps-aggregated results/RIBOSAM-norm-conds-by-reps-aggregated $RIBOSAM_PLUS_PDB_FASTA --nthreads=8
