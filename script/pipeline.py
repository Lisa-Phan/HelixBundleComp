"""
9/11/2025
Lisa P

Adapt original pipeline to run helix bundle design
for Casey on Quartz

ORIGINAL HOW TO RUN
1. Get the params needed for inpainting, conservation, and contigs
2. Update params for the steps
3. Run the first step, getting RFDiffusion output
4. Check output topolopy, pdbs that pass gets placed into new folder named MPNN_input
5. Run second step with MPNN. Check if the masking and indexing is correct, sometimes they are off
"""
########################
# IMPORTS
########################
import os
import sys
import textwrap
import json


########################
# PATHS CONFIGURATION
########################

GIT_DIR=r"/geode2/home/u020/cvanstap/HelixBundleComp"

CONFIG_DIR = r"/geode2/home/u020/cvanstap/rf_diffusion_all_atom/config"
SCRIPT_DIR = GIT_DIR+r"/script/utils" #accessory functions
SBATCH_MODE = False
LIGAND_MPNN_DIR = r"/geode2/home/u020/cvanstap/LigandMPNN"


########################
# TEMPLATE FILES
########################
#TODO: isolate alll relvant path variables to be changed and put it here
TEMPLATE_VARIABLE_DICT = {
    'email': 'cvanstap@iu.edu',
    'working_dir': GIT_DIR +  r"/Trial_run",
    'rfdiff': {'model': ''}
}

# INSTALL A SEPARATE CONDA ENVIRONMENT FOR BIOTITE
PYTHONS = {'biotite': '/N/u/cvanstap/Quartz/.conda/envs/biotite/bin/python'}

########################
# UTILITIES
########################
def generate_fixed_index_str(file):
    """
    generate conservation string to use in 
    LigandMPNN residue masking syntax

    Does not account for multiple chains jsonl
    return only last value
    """
    with open(file) as infile:
        json_obj = json.load(infile)
        for key in json_obj:
            chain_indices_dict = json_obj[key]
            for chain in chain_indices_dict:
                indices = json_obj[key][chain]
                print(f'keys {key}, chain {chain}, indices {indices}')
                conservation_str = ' '.join(f'{chain}{index}' for index in indices)
                print('conservation str: ', conservation_str)
    return conservation_str

#################################
# STEP 1 RFDiffusion
#################################

RFDiffusion_general_template = textwrap.dedent("""#!/bin/bash
#----------------------------------------------------------------------
#SBATCH -J rf_{job_name}
#SBATCH -o %x.out      # Name of stdout output file (%x expands to job name)
#SBATCH -e %x.err      # Name of stderr output file (%x expands to job name)
#SBATCH -p gpu-debug
#SBATCH -A r01589
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cvanstap@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=01:00:00
#----------------------------------------------------------------------
                                               
module load apptainer
which apptainer

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate SE3nv 

##############
# PATHS
##############
RFAA_PATH=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/
CONTAINER=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/rf_se3_diffusion.sif
INFERENCE=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/run_inference.py
CONFIG_FILE="{config_file}"
CONFIG_DIR=/geode2/home/u020/cvanstap/rf_diffusion_all_atom/config

HYDRA_FULL_ERROR=1

INPUT="{input_pdb}"
OUTPUT_DIR=RFAA_out/sample
mkdir -p $OUTPUT_DIR

#not sure why path can not be found, try changing directory to RFDiffusion
#as temporary fix
cd $RFAA_PATH


""")


RFDiffusion_run_command = textwrap.dedent("""
apptainer run --nv $CONTAINER -u $INFERENCE \\
--config-dir=$CONFIG_DIR \\
--config-name=$CONFIG_FILE \\
inference.deterministic=True \\
inference.output_prefix=$OUTPUT_DIR \\
inference.input_pdb=$INPUT \\
inference.ckpt_path=$MODEL \\
inference.num_designs=20 \\
inference.design_startnum=0 \\
""")

RFDiffusion_config = textwrap.dedent("""
defaults:
  - aa

diffuser:
  T: 200

inference:
  num_designs: 20
  model_runner: NRBStyleSelfCond
  ligand: 'Fe'

model:
  freeze_track_motif: True

contigmap:
  contigs: {contigs}
  inpaint_str: {inpaint_str}
  length: "{length}"

potentials:
  guiding_potentials: ["type:ligand_ncontacts,weight:1"] 
  guide_scale: 2
  guide_decay: cubic
""")

def step1_runRFDiffusion(config_file: str, 
                        job_name: str,  
                        inpaint_str: str,
                        contigs: str,
                        length: str,
                        input_pdb: str, 
                        runcommand: str = RFDiffusion_run_command,
                        config_dir = CONFIG_DIR
                        ):
    """
    Set up an RFDiffusion run job
    starting out with 
    1. structure 
    2. masking residue indices, or inpaint_str: parts where you want RFDiffusion to change
    3. contigs: parts where you want RFDiffusion to keep the same
    """
    #make config file 
    config_file_path = os.path.abspath(os.path.join(config_dir, f"inference/{config_file}"))
    print('config file path: ', config_file_path)
    with open(config_file_path, "w") as f:
        f.write(RFDiffusion_config.format(contigs=contigs, 
                                              inpaint_str=inpaint_str, 
                                              length=length))
        f.close()
    
    #make job file
    job_file = f"RFDiffusion_{job_name}.sh"
    with open(job_file, "w") as f:
        f.write(RFDiffusion_general_template.format(job_name=job_name,
                                                    email=TEMPLATE_VARIABLE_DICT['email'],
                                                    input_pdb=input_pdb, 
                                                    config_file = os.path.basename(config_file_path)))
        f.write(runcommand)
    #run the job
    if SBATCH_MODE:
      os.system(f"sbatch {job_file}")
    else:
      os.system(f"bash {job_file}")


#########################
# STEP 2 PROTEINMPNN
#########################       

ligandMPNN_template = textwrap.dedent("""#!/bin/bash
#----------------------------------------------------------------------
#SBATCH -J mpnn_{basename}  # Job name
#SBATCH -o %x.out      
#SBATCH -e %x.err
#SBATCH -p gpu-debug
#SBATCH -A r01589
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user={email}
#----------------------------------------------------------------------
                               
#initialize conda
eval "$(conda shell.bash hook)"
conda activate /N/u/cvanstap/Quartz/.conda/envs/ligandMPNN

LIGAND_MPNN_DIR=/geode2/home/u020/cvanstap/LigandMPNN

PDB_PATH={pdb}
BASENAME={basename}
JOB_NAME={job_name}
OUT_DIR={out_dir}/{basename}
                                                        
################
# Input
################

mkdir -p $OUT_DIR
""")

soluble_ligandMPNN_run_command = textwrap.dedent("""
python $LIGAND_MPNN_DIR/run.py \\
        --model_type "soluble_mpnn" \\
        --checkpoint_soluble_mpnn "$LIGAND_MPNN_DIR/model_params/solublempnn_v_48_010.pt" \\
        --pdb_path $PDB_PATH \\
        --out_folder $OUT_DIR \\
        --pack_side_chains 1 \\
        --number_of_packs_per_design 4 \\
        --pack_with_ligand_context 1 \\
        --fixed_residues "{fixed_residue_str}" \\
        --repack_everything 1 \\
        --batch_size 1 \\
        --number_of_batches 20
""")

def step2_runProteinMPNN(input_dir: str,
                        job_name: str,
                        conserved_res: str,
                        template_structure: str,
                        model_type: str = 'ligandMPNN_sol'):
    """
    First, make a jsonl file for position masking
    Implemented residue mapping is based on structure superposition
    This is automation, but manual check is recommended 
    Then run ProteinMPNN with fixed mappings

    By default run ligandMPNN_sol
    """
    os.makedirs(f'{job_name}_jsonl_dir', exist_ok=True)
    json_dir_realname = os.path.realpath(f'{job_name}_jsonl_dir')

    #call the script that makes the jsonl file
    os.system(f"{PYTHONS['the_other_biotite']} {SCRIPT_DIR}/create_fixed_residue_jsonl.py -i {input_dir} -o {job_name}_jsonl_dir -r {conserved_res} -f {template_structure}")

    runfile_dir = f'MPNN_runfile_dir'
    os.makedirs(runfile_dir, exist_ok=True)

    output_dir = f"MPNN_output"
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.realpath(output_dir)

    for file in os.listdir(input_dir):
        if file.endswith(".pdb"):
            job_name = os.path.basename(file).split(".")[0]
            pdb_path = os.path.join(input_dir, file)
            jsonl_path = os.path.join(json_dir_realname, f"{job_name}_fixed.jsonl")
            
            runfile = f"{runfile_dir}/{job_name}.sh"
            with open(runfile, "w") as f:

                if model_type == 'ligandMPNN_sol': 
                    fixed_residue_str = generate_fixed_index_str(jsonl_path)
                    #create a symlink to ligandMPNN params to current script directory
                    os.system(f'ln -sf {LIGAND_MPNN_DIR}/model_params ./')
                    f.write(ligandMPNN_template.format(pdb = pdb_path,
                                                        job_name = job_name,
                                                        email = TEMPLATE_VARIABLE_DICT['email'], 
                                                        out_dir = output_dir, 
                                                        basename = job_name))
                    f.write(soluble_ligandMPNN_run_command.format(fixed_residue_str = fixed_residue_str))
                
                else: 
                    print('model_type {model_type} is not supported. Exiting')
                    exit()


    #run the jobs
    for file in os.listdir(runfile_dir):
        if file.endswith(".sh"):
            if SBATCH_MODE: 
              os.system(f"sbatch {runfile_dir}/{file}")
            else:
              os.system(f"bash {runfile_dir}/{file}")
            #wait
            os.system("sleep 1")

########################
# TODO: update ended here on 9/11/25
# Update remaining component tomorrow
#########################


########################################
# STEP3 run Alphafold
########################################

ALPHAFOLD3_TEMPLATE = textwrap.dedent("""#!/bin/bash
#----------------------------------------------------------------------
#SBATCH -o %x.out         # Name of stdout output file
#SBATCH -e %x.err         # Name of stderr error file
#SBATCH -p gpu-a100                # Queue (partition) name
#SBATCH -N 1                       # Total # of nodes
#SBATCH -t 2:00:00                # Run time (hh:mm:ss)
#SBATCH -A CHE23010              # Allocation name
#----------------------------------------------------------------------

# Load required modules
module unload xalt
module use /scratch/tacc/apps/bio/alphafold3/modulefiles
module load alphafold3/3.0.1-ctr.lua
module load tacc-apptainer
echo loaded modules
module list

# Set environment variable definitions to point to your input, output, and model parameters directories:
export AF3_INPUT_DIR={INPUT_DIR}
export AF3_OUTPUT_DIR={OUTPUT_DIR}
export AF3_MODEL_PARAMETERS_DIR=/work/09069/dhp563/ls6/alphafold3_weight/

# Run AlphaFold3 container 
apptainer exec --nv -B $AF3_INPUT_DIR:/root/af_input -B $AF3_OUTPUT_DIR:/root/af_output -B $AF3_MODEL_PARAMETERS_DIR:/root/models -B $AF3_DATABASES_DIR:/root/public_databases $AF3_IMAGE python $AF3_CODE_DIR/run_alphafold.py --json_path=/root/af_input/input.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output
""")

GENERIC_SUBMIT_TEMPLATE = textwrap.dedent("""#!/bin/bash
#SBATCH -A CHE23010
#SBATCH -p gpu-a100     # Queue name
#SBATCH -t 00:10:00
#SBATCH -N 1

""")

#Check this
def step1():
    inpaint_str = '["197,200,205,209,213,216,217,220,221,224,225,228,236,239,243,245,246,247,249,250,253,254,271,272,274,275,279,282,285,286,288,289,292,293,295,296,297,304,305,306,308,309,310,312,313,316,317,319,320,321,324,325,328,333,334,342,343,345,346,349,350,353,357,358,360,361,364,365,368"]'
    contigs = '["B196-231,2-4,B234-260,5-10,B267-297,7-10,B304-368"]'
    length = '169-179'
    input_pdb = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/input/NOR_chainB_heme_copper.pdb"
    
    print('parameters: ')
    print('input: ', input_pdb)
    
    step1_runRFDiffusion("cNOR.yaml", 
                         "cNOR", 
                         inpaint_str, 
                         contigs, 
                         length, 
                         input_pdb)


#conserved resi is a fairly wide selection from a 
#modified 75% conservation threshold, note that this is a
#much larger selection compared to the original HCO selection
def step2():
    input_dir = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/MPNN_input"
    job_name = "cNOR"
    conserved_res = "198,202,203,206,207,209,210,211,212,214,215,252,255,256,257,258,259,260,262,263,264,266,268,269,273,276,277,280,283,311,322,323,326,327,330,332,336,338,339,340,341,342,343,346,347,348,352,355,356"
    step2_runProteinMPNN(input_dir, 
                         job_name, 
                         conserved_res, 
                         model_type='ligandMPNN_sol',
                         template_structure=f"{TEMPLATE_VARIABLE_DICT['working_dir']}/input/NOR_chainB_heme_copper.pdb")

def step3():
    """
    Make AF3 json from MPNN output
    """
    #nested directory of MPNN_output/jobname_rfdiffusion_number/seqs/fastafile
    input_dir = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/MPNN_output"
    os.makedirs(f"{TEMPLATE_VARIABLE_DICT['working_dir']}/AF3_json", exist_ok=True)
    #get all files under nested directory ending with .fa
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fa"):
                #get the full path of the file
                file_path = os.path.join(root, file)
                #run script to make AF3 json
                #single option is for if you want to submit one file on the server
                # multi is for when you want to locally on TACC
                os.system(f'cd AF3_json && bash {SCRIPT_DIR}/make_af3json.sh {file_path} multi')
    
    #create alphafold3 directories to run job
    os.makedirs(f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_input", exist_ok=True)
    os.makedirs(f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_output", exist_ok=True)
    os.makedirs(f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_runfiles", exist_ok=True)

    for file in os.listdir(f"{TEMPLATE_VARIABLE_DICT['working_dir']}/AF3_json"):
        if file.endswith(".json"):
            #make a directory under alphafold3_input with the same name as the json file
            basename = os.path.splitext(file)[0]
            
            input_dir = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_input/{basename}"
            output_dir = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_output/{basename}"

            os.makedirs(input_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)
            #move the json file to the input directory
            os.rename(f"{TEMPLATE_VARIABLE_DICT['working_dir']}/AF3_json/{file}", f"{input_dir}/input.json")

            #make a run file for alphafold3
            runfile = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_runfiles/{basename}.sh"
            with open(runfile, "w") as f:
                f.write(ALPHAFOLD3_TEMPLATE.format(INPUT_DIR=input_dir, 
                                                   OUTPUT_DIR=output_dir))

    #create a sbatch file that submits the alphafold3 jobs
    #break into chunks of 30 jobs at a time
    runfiles_dir = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_runfiles"
    file_list = [f for f in os.listdir(runfiles_dir) if f.endswith('.sh') and not f.startswith('submit_')]
    for i in range(0, len(file_list), 30):
        chunk = file_list[i:i+30]
        sbatch_file = f"{TEMPLATE_VARIABLE_DICT['working_dir']}/alphafold3_runfiles/submit_chunk_{i//30}.sh"
        with open(sbatch_file, "w") as f:
            f.write(GENERIC_SUBMIT_TEMPLATE.format(chunk_name=f"chunk_{i//30}"))
            for file in chunk:
                f.write(f"sbatch -N 1 -t 2:00:00 {runfiles_dir}/{file}\n")
            f.close()
                

if __name__ == "__main__":
    step2()
