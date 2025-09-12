#!/bin/bash
#SBATCH -J mpnn  # Job name
#SBATCH -o %x.out      # Name of stdout output file (%j expands to job ID)
#SBATCH -e %x.err      # Name of stderr output file (%j expands to job ID)
#SBATCH -A CHE23010
#SBATCH -p gpu-a100     # Queue name
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=dhp563@utexas.edu

eval "$(conda shell.bash hook)"
conda activate /work/09069/dhp563/miniconda_installation_20230606/envs/ligandmpnn_env #change this

LIGAND_MPNN_DIR=/work/09069/dhp563/ls6/LigandMPNN

PDB_PATH=/work/09069/dhp563/ls6/HelixBundleComp/structure_query/input/2INP_core_diiron.pdb
BASENAME=2INP_core_simple
JOB_NAME=2INP_core_simple
OUT_DIR=./sample/2INP_core_withactive

#create symlink for model param because program cannot recognize provided path

ln -s /work/09069/dhp563/ls6/LigandMPNN/model_params ./ 

################
# Input
################

mkdir -p $OUT_DIR




python $LIGAND_MPNN_DIR/run.py \
        --model_type "soluble_mpnn" \
        --checkpoint_soluble_mpnn $LIGAND_MPNN_DIR/model_params/solublempnn_v_48_010.pt \
        --pdb_path $PDB_PATH \
        --fixed_residues "B108 B138 B141 B199 B233 B236" \
        --out_folder $OUT_DIR \
        --pack_side_chains 1 \
        --number_of_packs_per_design 4 \
        --repack_everything 1 \
        --batch_size 1 \
        --number_of_batches 20 \
        --omit_AA='MC'
