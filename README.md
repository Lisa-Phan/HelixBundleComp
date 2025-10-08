
7/18/2025

# Demo comp pipeline to query and design helix bundle protein

## Step 1: Get candidate structures


Here we query a few representatives PDB for diiron core helix bundle, based on Casey's suggestion. One for sMMO, one for T4MO, one for RNR, one for PHPi
Search will be done on Foldseek's PDB and AFUni30 database, with CA-backbone template and TSV format


## Step 2: Filter search


Filter hit criteria to ensure that reported structures are likely homologs


## Step 3: Get core sequence


Run foldseek's result2msa to get core the aligned sequence that define hit structure's helix bundle. Combine this 
with the filtered result to provide inputs for the next step.



