#!/bin/bash

#SBATCH -J test_mpnn
#SBATCH -p gpu-debug
#SBATCH -A r01589
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cvanstap@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:30:00

module load conda
conda activate proteinMPNN

PROTEIN_MPNN_PATH=/geode2/home/u020/cvanstap/ProteinMPNN
OUTFOLDER=./test_proteinMPNN

#Create sym link to model params in current dir
ln -s $PROTEIN_MPNN_PATH/model_params ./
mkdir -p $OUTFOLDER

#test the default example
python $PROTEIN_MPNN_PATH/run.py \
        --seed 111 \
        --pdb_path $PROTEIN_MPNN_PATH/inputs/1BC8.pdb \
        --out_folder $OUTFOLDER

rm ./model_params
                       
