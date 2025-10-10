#!/bin/bash

#using mmseqs to find similar sequences to a starting point
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate /stor/work/YiLu/conda/miniconda3/envs/mmseqs

mkdir database

mmseqs databases UniProtKB/Swiss-Prot database/UniProt databasetmp
mmseqs createdb input/* database/inputdb
mmseqs createindex database/inputdb indextmp

#create the tsv with easy-search, optional
mmseqs easy-search input/* database/UniProt output/UniProt_hit.tsv searchtmp \
-s 7.5 \
--format-mode 4 \
--format-output 'query,target,evalue,pident,qstart,qend,tstart,tend,qaln,taln'

#create an alignment
mmseqs search database/inputdb database/UniProt database/UniProt_hit searchtmp -a -s 7.5
mmseqs result2msa database/inputdb database/UniProt database/UniProt_hit output/UniProt_hit_msa.fasta --msa-format-mode 4
