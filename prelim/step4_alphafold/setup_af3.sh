#!/bin/bash 

# Lisa P
# 11/2/25
# run_alphafold.sh
# 
# main script to control directory setup to run alphafold

#################
#SLURM TEMPLATE
#################
fill_slurm_template () {
    AF_INDIR=$(realpath $1)
    AF_OUTDIR=$(realpath $2)
    AF_JOBNAME=$3

    cat <<EOF
#!/bin/bash

#SBATCH -J $AF_JOBNAME
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -A r01589
#SBATCH -p hopper
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --gpus-per-node 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cvanstap@iu.edu
#SBATCH --time=1:00:00

module load alphafold/3.0.0

export AF_MODELS_DIR=/N/slate/$USER/alphafold_models
export AF_INPUTDIR=$AF_INDIR
export AF_OUTDIR=$AF_OUTDIR

#cp alphafold_input.json $AF_INPDIR
cp /N/soft/rhel8/alphafold/3.0.0/example/run_alphafold.py $AF_INDIR

apptainer exec \\
    --nv \\
    --bind $AF_INPUTDIR:/root/af_input \\
    --bind $AF_OUTDIR:/root/af_output \\
    --bind $AF_MODELS_DIR:/root/models \\
    $AF_CONTAINER \\
    python /root/af_input/run_alphafold.py \\
    --json_path=/root/af_input/alphafold_input.json \\
    --model_dir=/root/models \\
    --db_dir=/root/public_databases \\
    --output_dir=/root/af_output
EOF
}




####################
# Step 1, make jsons
####################
WORKINGDIR=$(pwd)
LOG=output/log.txt
touch $LOG

mkdir tmp
cd tmp
bash $WORKINGDIR/make_af3json.sh $WORKINGDIR/input/*.fa multi >> ../output/log.txt

#####################
# Step2, set up SLURM runfiles
#####################
cd $WORKINGDIR
mkdir -p {slurm_runfiles,alphafold3_input,alphafold3_output}

#for each .json file, create a subdir under _input and _output
for file in $(ls tmp/*.json); do 
    #get basenname
    fileBaseName=$(basename $file)
    
    #get name without .json
    jobName=${fileBaseName%.*}

    #make subdirectories
    inDir=$(realpath alphafold3_input/$jobName)
    outDir=$(realpath alphafold3_output/$jobName)

    mkdir $inDir
    mkdir $outDir
    
    #create slurm template
    slurmFile=slurm_runfiles/submit_$jobName.sh
    echo Making slurm file $slurmFile >> $LOG
    fill_slurm_template "$inDir" "$outDir" "$jobName" >> $slurmFile

    #move json into input.json
    mv $(realpath $file) $inDir/alphafold_input.json

done
