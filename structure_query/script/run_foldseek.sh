#!/bin/bash

# Foldseek download common databases
conda activate foldseek

#this will take a few minutes
foldseek databases Alphafold/UniProt50-minimal AFUniprot50min tmp 
foldseek databases PDB PDB_db PDBtmp

#make inputDB
foldseek createdb input inputDB

#do a format mode 4 for table specs output
# and format mode 5 search for backbone
#since we are doing topology sampling, there is 
#no need to be exhaustive yet
foldseek easy-search inputDB PDB_db PDB_search_result search_tmp \
--cluster-search 0 \
--format-mode 4 \
--format-output query,target,fident,alnlen,mismatch,evalue,bits,tseq,lddt,qtmscore,ttmscore,rmsd,prob

foldseek easy-search database/inputDB database/PDB_db PDB_search_format5_ search_tmp --format-mode 5

foldseek easy-search inputDB AFUniprot50min AFUni50min_search_result AFUni_search_tmp \
--cluster-search 0 \
--format-mode 4 \
--format-output query,target,fident,alnlen,mismatch,evalue,bits,tseq,lddt,qtmscore,ttmscore,rmsd,prob

foldseek easy-search database/inputDB database/AFUniprot50min AFUni50min_format5_ tmp --format-mode 5
