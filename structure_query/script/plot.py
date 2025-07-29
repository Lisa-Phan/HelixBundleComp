"""
7/22/2025
Lisa P. 

Simple data analysis to see what the RMSD 
hit distribution is
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns   

FOLDSEEK_HIT = r"/stor/scratch/YiLu/dhp563/HelixBundleComp/structure_query/AFUni50min_search_result.tsv"

def plot_data(file):
    df = pd.read_csv(file, sep="\t")
    print(df.head())
    #add an RMSD filter of 2
    df = df[df['rmsd'] < 2]
    #plot RMSD distribution among queries as histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(data = df, x='rmsd', hue='query')
    plt.title("RMSD Distribution of <2 Angstrom hits from \nAlphaFold Uniprot 50 hits")
    plt.tight_layout()
    plt.savefig("rmsd_distribution.png")

plot_data(FOLDSEEK_HIT)
