#!/bin/bash

# Make input template for alphafold3
# Take MPNN output and place into alphafold3 json template
# Usage: ./make_af3json.sh <MPNN_FASTA> <single|multi>

fill_template () {
    # Print JSON with proper quoting
    JOB=$1
    SEQUENCE=$2

    cat <<EOF
    {
    "name": "$JOB",
    "modelSeeds": [],
    "sequences": [
        {
        "proteinChain": {
            "sequence": "$SEQUENCE",
            "count": 1
        }
        },
        {
        "ion": {
            "ion": "FE",
            "count": 2
        }
        }
    ],
    "dialect": "alphafoldserver",
    "version": 1
    }
EOF
}


singlefile_template_function() {
    FASTA_PATH=$1
    #clean if file exists
    if [ -f out.json ]; then
        rm out.json
    fi

    # Initialize variables
    JOB_NAME=""
    SEQUENCE=""

    echo "[" >> out.json
    while read -r line; do
        if [[ $line == ">"* ]]; then
            JOB_NAME="${line#>}"
            #use the id and the starting name as job name
            JOB_NAME=$(echo $JOB_NAME | cut -d',' -f1,2 | sed 's/, id=/_id_/g')
            echo job_name: $JOB_NAME
        else
            SEQUENCE="${line// /}"
            echo sequence: $SEQUENCE
            fill_template $JOB_NAME $SEQUENCE >> out.json
            # Add a comma if it's not the last record
            echo "," >> out.json        
        fi
    done < "$FASTA_PATH"
    # Remove the last comma
    sed -i '$ s/,$//' out.json
    echo "]" >> out.json
}

multifile_template_function() {
    FASTA_PATH=$1
    #clean if file exists
    if [ -f out.json ]; then
        rm out.json
    fi

    # Initialize variables
    JOB_NAME=""
    SEQUENCE=""
    OUTFILE_NAME="default.json"


    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            JOB_NAME="${line#>}"
            JOB_NAME=$(echo "$JOB_NAME" | cut -d',' -f1,2 | sed 's/, id=/_id_/g')
            echo "job_name: $JOB_NAME"
            OUTFILE_NAME="${JOB_NAME}.json"
            echo "[" >> "$OUTFILE_NAME"
        else
            if [[ -z "$OUTFILE_NAME" ]]; then
                echo "Error: OUTFILE_NAME not set before writing. Line: $line" >&2
                continue
            fi
            SEQUENCE="${line// /}"
            echo "sequence: $SEQUENCE"
            fill_template "$JOB_NAME" "$SEQUENCE" >> "$OUTFILE_NAME"
            echo "]" >> "$OUTFILE_NAME"
        fi
    done < <(tail -n +3 "$FASTA_PATH")
}




FASTA=$1
ARG=$2

if [ $ARG == 'multi' ]; then
    multifile_template_function $FASTA
elif [ $ARG == 'single' ]; then
    singlefile_template_function $FASTA
else
    echo 'MPNN to AF3json template generator'
    echo 'Usage: make_af3json <fasta> <single/multi>'
fi

