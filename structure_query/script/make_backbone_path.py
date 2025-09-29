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
AFUni_filtered = r"./prelim_filtered_result/AFUniprot_filtered_result.csv"
PDB_filtered = r"./prelim_filtered_result/PDB_filtered_result.csv"

AFUni_backbone_path = r"../hit_back_bone/AFUni_hit_backbone"
PDB_backbone_path = r"../hit_back_bone/PDB_hit_backbone"

################
# FUNCTION 
################
def path_mapping(hit_pdb_string: str, template_string: str, file_list: list, dir: str):
    """
    return the fullpath of the pdb mentioned in substring, if they are included
    somewhere in filelist
    """
    ans = []
    for file in file_list: 
        #if the hit pdb is a substring
        if hit_pdb_string in file:
            if template_string in file: 
            #check for a template match as well
                ans.append(os.path.join(dir, file))
    
    return ans[0]

QUERY_PDB_MAP = {
    "1W69_core" : "../input/1W69_core.pdb",
    "2INP_core" : "../input/2INP_core.pdb",
    "5TDU_core" : "../input/5TDU_core.pdb",
    "6YD0_core" : "../input/6YD0_core.pdb"
}

###############
# MAIN
###############

AFUni_filelist = os.listdir(AFUni_backbone_path)
PDB_filelist = os.listdir(PDB_backbone_path)

AFUni_table = pd.read_csv(AFUni_filtered)
AFUni_table['hit_path'] = AFUni_table.apply(lambda row: path_mapping(row['target'], 
                                                                     row['query'], 
                                                                     AFUni_filelist, 
                                                                     AFUni_backbone_path),
                                                                     axis=1)

PDB_table = pd.read_csv(PDB_filtered)
PDB_table['hit_path'] = PDB_table.apply(lambda row: path_mapping(row['target'],
                                                                 row['query'],
                                                                 PDB_filelist, 
                                                                 PDB_backbone_path), 
                                                                 axis=1)

#first two columns are query target
joined = pd.concat([AFUni_table, PDB_table])
joined['query_path'] = joined['query'].apply(lambda x: QUERY_PDB_MAP[x])

joined[['query_path', 'hit_path']].to_csv('Path_pairing_for_hit_trimming.csv', index=False)




