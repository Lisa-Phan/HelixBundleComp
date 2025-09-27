"""
9/26/25
Lisa P

Create a path file for backbone filtering
TODO: discuss what the outcome should be when 
a hit matches with multiple queries
"""
################
# IMPORTS
################
import pandas as pd
import numpy as np
import os

################
# PATHS
################
#prefiltered by fident <= 0.7 and prob > 0.3 for PDB and 0.5 for AFUniprot 
AFUni_filtered = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/script/prelim_filtered_result/AFUniprot_filtered_result.csv"
PDB_filtered = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/script/prelim_filtered_result/PDB_filtered_result.csv"

AFUni_backbone_path = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/hit_back_bone/AFUni_hit_backbone"
PDB_backbone_path = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/hit_back_bone/PDB_hit_backbone"

################
# FUNCTION 
################
def path_mapping(sub_string: str, file_list: list):
    """
    return the fullpath of the pdb mentioned in substring, if they are included
    somewhere in filelist
    """
    ans = []
    for file in file_list: 
        if sub_string in file: 
            ans.append(file)
    
    print(ans)


###############
# MAIN
###############

AFUni_filelist = os.listdir(AFUni_backbone_path)
PDB_filelist = os.listdir(PDB_backbone_path)

AFUni_table = pd.read_csv(AFUni_filtered)
AFUni_table['path'] = AFUni_table['target'].apply(lambda x: path_mapping(x, AFUni_filelist))

PDB_table = pd.read_csv(PDB_filtered)
PDB_table['path'] = PDB_table['target'].apply(lambda x: path_mapping(x, PDB_filelist))

#first two columns are query target
joined = pd.concat([AFUni_table, PDB_table])
print(joined.head())



