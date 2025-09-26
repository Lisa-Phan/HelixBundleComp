"""
9/13/2025
Lisa P
Foldseek filtering script 

"""

import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

foldseek_file = r"PDB_search_result.tsv"

result = pd.read_csv(foldseek_file, sep='\t')

result = result.astype({'fident': 'float64',
                        'lddt' : 'float64', 
                        'rmsd': 'float64', 
                        'prob': 'float64'})

#plot histogram of probability 
#prob_hist = sns.histplot(result, x='prob')
#plt.savefig('prob_dist.png')


fident_hist = sns.histplot(result, x='fident')
plt.savefig('fident_dist.png')

#filter numbers
filtered = result[(result['prob'] <= 0.8) & (result['fident'] >= 70)]

