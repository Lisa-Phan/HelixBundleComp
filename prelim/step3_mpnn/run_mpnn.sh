#!/bin/bash
#SBATCH -J mpnn  # Job name
#SBATCH -o ./output/%x.out      # Name of stdout output file (%j expands to job ID)
#SBATCH -e ./output/%x.err      # Name of stderr output file (%j expands to job ID)
#SBATCH -A CHE23010
#SBATCH -p gpu-a100     # Queue name
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=dhp563@utexas.edu

eval "$(conda shell.bash hook)"
conda activate /work/09069/dhp563/miniconda_installation_20230606/envs/ligandmpnn_env #change this

LIGAND_MPNN_DIR=/work/09069/dhp563/ls6/LigandMPNN

PDB_PATH=./input/1W69_core_diiron.pdb
BASENAME=RNR_core
JOB_NAME=RNR_core_test1
OUT_DIR=output/1W69_core

################
# Input
################

mkdir -p $OUT_DIR

python $LIGAND_MPNN_DIR/run.py \
        --model_type "soluble_mpnn" \
        --checkpoint_soluble_mpnn "$LIGAND_MPNN_DIR/model_params/solublempnn_v_48_010.pt" \
        --pdb_path $PDB_PATH \
        --out_folder $OUT_DIR \
        --pack_side_chains 1 \
        --number_of_packs_per_design 4 \
        --pack_with_ligand_context 1 \
        --fixed_residues "A143 A135 A171 A128 A146 A167 A174 A140 A178 A121 A164 A125 A134 A198 A215 A147 A133 A175 A242 A136 A138 A159 A176 A209 A179 A188 A238 A256 A257 A137 A144 A145 A151 A196 A269 A142 A156 A172 A236 A245 A246 A163 A226 A234 A235 A275 A129 A199 A208 A211 A231 A249 A261 A264 A266 A267 A274 A169 A173 A195 A230 A233 A248 A262 A268 A271 A276" \
        --repack_everything 1 \
        --batch_size 1 \
        --number_of_batches 20 \
        --omit_AA='MC'
