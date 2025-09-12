#!/bin/bash
#SBATCH -J rf_5B1B_extended_inpaint_aggressive  # Job name
#SBATCH -o %x.out      # Name of stdout output file (%j expands to job ID)
#SBATCH -e %x.err      # Name of stderr output file (%j expands to job ID)
#SBATCH -A CHE23010
#SBATCH -p gpu-a100-dev     # Queue name, needs to be gpu
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=dhp563@my.utexas.edu

module load tacc-apptainer
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate SE3nv 
module load tacc-apptainer

OUTPUT_PREFIX=RFDiffusion_output/5B1B_extended_inpaint_aggressive
INPUT=/work/09069/dhp563/ls6/HCO_modelling/branch_pipeline/input/5B1B_chainA_heme_copper_cutheme.pdb
#MODELS="/work/09069/dhp563/ls6/RFdiffusion/models"

#########################
#for all atom RFDiffusion
#########################
ALL_ATOM_SIF="/work/09069/dhp563/ls6/rf_diffusion_all_atom/rf_se3_diffusion.sif"
ALL_ATOM_INFERENCE="/work/09069/dhp563/ls6/rf_diffusion_all_atom/run_inference.py"
MODEL="/work/09069/dhp563/ls6/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt"
CONFIG_FILE="5B1B_extended_inpaint_aggressive.yaml"
CONFIG_DIR="/work/09069/dhp563/ls6/rf_diffusion_all_atom/config"

#for stack tracing
HYDRA_FULL_ERROR=1

module load cuda

apptainer run --nv $ALL_ATOM_SIF -u $ALL_ATOM_INFERENCE \
--config-dir=$CONFIG_DIR \
--config-name=$CONFIG_FILE \
inference.deterministic=True \
inference.output_prefix=$OUTPUT_PREFIX \
inference.input_pdb=$INPUT \
inference.ckpt_path=$MODEL \
inference.num_designs=20 \
inference.design_startnum=0 \
