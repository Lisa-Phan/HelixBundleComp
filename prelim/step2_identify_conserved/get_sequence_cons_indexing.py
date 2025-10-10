"""
10/10/2025

Lisa P. 
Re-write script to pick out residues to mask and to conserve
Script originally run on rentcomp

Steps: 
1. Generate an msa from starting fasta
2. Use script to create conservation map on corresponding pdb and find things to conserve
"""
import os
import sys
import textwrap
import biotite
import matplotlib.pyplot as plt
from biotite.sequence.io.fasta import FastaFile
from biotite.sequence.align import Alignment
import numpy as np


BIN = {'mmseqs': '/stor/work/YiLu/conda/miniconda3/envs/mmseqs/bin/mmseqs'}

###########################
# For helix bundle comp
###########################
#input fasta is the same as pdb sequence
#both generated from the truncated structure
INPUT_FA = r'/stor/scratch/YiLu/dhp563/HelixBundleComp/prelim/step1_query_sequence/input/1W69_core.fasta'

#HMM_FASTA = r'/stor/scratch/YiLu/dhp563/HCO_search/literature/HCO/Sequence_database/HCO-family-sequences_unique/CNOR_sequences_unique.txt'
ALIGNMENT = r'/stor/scratch/YiLu/dhp563/HelixBundleComp/prelim/step1_query_sequence/output/UniProt_hit_msa.fasta'

#indices making up the helix around the heme-copper binding
#COMPLETE_INDICES = np.arange(195, 368)

#conservation threshold for residues to be marked as conserved
CONS_THRESHOLD = 0.6


#residue number that started the protein sequence. = 0 if starting from res 1
#for rnr bundle this is 120
OFFSET = 120

#a little unnecessary compared to directly building a msa, but 
#this should still work
MSA_gen_template = """
#!/bin/bash
conda activate mmseqs

QUERY={query}
HMM_FASTA={hmm_fasta}
DATABASE="/stor/scratch/YiLu/dhp563/HCO_search/mmseqs_search_hco/databases"

mmseqs createdb $QUERY $DATABASE/{query_name}DB
mmseqs createdb $HMM_FASTA $DATABASE/{hmm_name}DB 
mmseqs search $DATABASE/{query_name}DB $DATABASE/{hmm_name}DB {hmm_name}_resultDB $DATABASE/tmp --num-iterations 4 --exhaustive-search --max-seq-id 0.5 
mmseqs result2msa $DATABASE/{query_name}DB $DATABASE/{hmm_name}DB {hmm_name}_resultDB {hmm_name}_msa
mmseqs unpackdb {hmm_name}_msa {hmm_name}_msa_unpacked_outdir
"""

def generate_msa(query, hmm_fasta, query_name, hmm_name):
    """
    Generate MSA using a starting PDB and sequences listed in a given HMM file
    """
    msa_script = MSA_gen_template.format(query=query, 
                                        hmm_fasta=hmm_fasta,
                                        query_name=query_name, 
                                        hmm_name=hmm_name)
    msa_script_path = f"{query_name}_msa.sh"
    with open(msa_script_path, 'w') as f:
        f.write(textwrap.dedent(msa_script))
    
    #does not work, need to run from command line
    #os.system(f'bash {msa_script_path}')

    return f'{hmm_name}_msa'


def get_index_dictionary(msa_file, indices=[], refname=None):
    """
    Read stacked fasta file in msa format and 
    get frequency for counts of residues mentioned at indices
    

    """
    aligned_seqs = FastaFile().read(msa_file)
    sequence_list = aligned_seqs.lines[1::2]
    sequence_header = aligned_seqs.lines[0::2]
    seq_dict = dict(zip(sequence_header, sequence_list))
    
    index_dict = {}
    #first one is the reference
    if refname and len(indices) == 0 : 
        ref_seq = aligned_seqs[refname]
        indices = range(0, len(ref_seq))

    for index in indices: 
        for key in seq_dict:
            try: 
                residue = seq_dict[key][index - 1]
                try:
                    index_dict[index].append(residue)
                except KeyError:
                    index_dict[index] = [residue]
            except IndexError:
                print(f'IndexError for index {index} and key {key}')
                continue
    return index_dict

def compute_frequency(frequency_dict, count_blank=False):
    """
    Compute frequency of residues in frequency dictionary
    return dictionary of indices : {residue : frequency}
    """
    for index_key in frequency_dict:
        all_residues = list(set(frequency_dict[index_key]))
        residue_count = {}
        for residue in all_residues:
            count = frequency_dict[index_key].count(residue)
            if count_blank:
                frequency = count / len(frequency_dict[index_key])
            else:
                non_blank = len(frequency_dict[index_key]) - frequency_dict[index_key].count('-')
                frequency = count / non_blank
            residue_count[residue] = frequency
        frequency_dict[index_key] = residue_count
    return frequency_dict

def calculate_conservation(frequency_dict, reference_fasta):
    """
    Calculate conservation based on reference for each residue
    """
    conservation = {}
    ref_seq = FastaFile().read(reference_fasta).lines[1]
    for index in frequency_dict:
        ref_residue = ref_seq[index - 1]
        #print('ref residue:', ref_residue)  
        try:
            conservation[index] = frequency_dict[index][ref_residue]
        except KeyError:
            conservation[index] = 0
    return conservation

def plot_conservation_distribution(input_dict):
    """
    Plot the conservation distribution 
    """
    plot_values = list(input_dict.values())
    plt.hist(plot_values, bins=30)
    plt.savefig('conservation_dist.png')


def main():

    msa_file = generate_msa(
        query=INPUT_FA,
        hmm_fasta=HMM_FASTA,
        query_name='input_cNOR',
        hmm_name='cNOR_hmm')
    
    index_dict = get_index_dictionary(msa_file, COMPLETE_INDICES)
    frequency_dict = compute_frequency(index_dict, count_blank=True)
    conservation = calculate_conservation(frequency_dict, INPUT_FA)
    
    print('Conservation:', conservation)

    sorted_conservation = dict(sorted(conservation.items(), key=lambda item: item[1], reverse=True))
    filtered_conservation = {k: v for k, v in sorted_conservation.items() if v > CONS_THRESHOLD}
    sorted_keys = list(filtered_conservation.keys())
    sorted_keys.sort()
    
    #print in pymol residue selection syntax
    selection_str = '+'.join(str(index) for index in sorted_keys)
    print(selection_str)

def main2(): 

    index_dict = get_index_dictionary(msa_file = ALIGNMENT,
                                      refname = '1W69_core')
    frequency_dict = compute_frequency(index_dict, count_blank=True)
    conservation = calculate_conservation(frequency_dict, INPUT_FA)
    
    #plot_conservation_distribution(conservation)
    #distribution shows something that is mostly normal-like
    #based on graph, pick 0.6 as threshold

    #TODO: add residue conservation percentage in bfactor column

    sorted_conservation = dict(sorted(conservation.items(), key=lambda item: item[1], reverse=True))
    filtered_conservation = {k: v for k, v in sorted_conservation.items() if v > CONS_THRESHOLD}
    sorted_keys = list(filtered_conservation.keys())

    selection_str = '+'.join(str(index + OFFSET) for index in sorted_keys)
    print(selection_str)


if __name__ == '__main__':
    main2()

