#!/bin/bash

#SBATCH -J test_rfaa
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

module load apptainer
which apptainer

RFAA_PATH=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/
CONTAINER=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/rf_se3_diffusion.sif
INFERENCE=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/run_inference.py


OUTPUT_DIR=RFAA_out/sample
mkdir -p $OUTPUT_DIR

#not sure why path can not be found, try changing directory to RFDiffusion
#as temporary fix
cd $RFAA_PATH

apptainer run --nv $CONTAINER -u $INFERENCE \
inference.deterministic=True \
diffuser.T=100 \
inference.output_prefix=$OUTPUT_DIR \
inference.input_pdb=$RFAA_PATH/input/7v11.pdb \
contigmap.contigs=[\'150-150\'] \
inference.ligand=OQO \
inference.num_designs=1 \
inference.design_startnum=0
