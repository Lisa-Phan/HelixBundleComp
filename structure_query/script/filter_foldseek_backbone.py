"""
Lisa P. 
9/26/25

Foldseek structure filtering to get 
the overlapping region for hits

TODO:
determine template structure
read in both backbone from template and hit 
identify matching positions of start and end on the hit
return truncated sequence and name for hit 
Result would be a TSV with the following column 

<name_of_hit> <name_of_query> <truncated_sequence> <rmsd_using_truncated_sequence>

Also write a set of truncated backbone in separate directory

#NOTE: sequence alignment would not work in this case
because these are distantly related structure homolog
#A better way to check would be secondary structure string

#This works in principle but not with the backbone position provided by 
foldseek. 

If running ProteinMPNN, then the complete backbone is not needed
Likely best if issue can be simplified via just coordinate matching for now.
"""
import argparse
import biotite
import os
import numpy as np

import biotite.sequence as seq
import biotite.sequence.align as align
from biotite.structure.io.pdb import PDBFile
from biotite.structure import distance
from biotite.structure import Atom, AtomArray

import argparse

########################
# UNUSED
########################

def sequence_alignment(str_query: str, str_template: str) -> None: 
    """
    Read in sequences for trial alignment
    Serve as the additional check in addition to 
    brute force distance indexing check
    """
    protein_query = seq.ProteinSequence(str_query)
    protein_templ = seq.ProteinSequence(str_template)
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    alignment = align.align_optimal(protein_query, 
                                    protein_templ, 
                                    matrix, 
                                    gap_penalty = -8)
    
    for alignment_char in alignment: 
        print(alignment_char)

def get_pydssp_annotation(structure: AtomArray) -> str:
    """
    Get structure annotation through pydssp
    pydssp expects a PDB with only C, O, N, CA, so temporarily write a file for that

    """
    #get C, O, N, Ca
    backbone = structure[np.isin(structure.atom_name, ['CA', 'N', 'C', 'O'])]
    
    if os.path.exists('tmp.pdb'):
        os.system('rm tmp.pdb')

    if os.path.exists('dssp_anno.txt'):
        os.system('rm dssp_anno.txt')

    #write a temporary structure to be used
    #TODO: find a less brute way to compute this
    tmpOutFile = PDBFile()
    tmpOutFile.set_structure(backbone)
    tmpOutFile.write('tmp.pdb')

    #read this with pydssp
    os.system(f'pydssp tmp.pdb -o dssp_anno.txt')
    
    #get string annotation
    with open('dssp_anno.txt', 'r') as inFile:
        for line in inFile: 
            annotation = line.split()[0]
            return annotation

def get_helix_index(scs_annotation: str, type: str = 'start') -> int:
    """
    return residue index for either start of first helix or
    end of last helix
    """
    if type == 'start':
        for index, letter in enumerate([char for char in scs_annotation]):
            if letter == 'H':
            #found first index, return residue number, which starts from 1
                return index+1
            
    elif type == 'end':
        index = -1
        tail = scs_annotation[index]
        while tail != 'H':
            index -= 1
            tail = scs_annotation[index]
        return len(scs_annotation) - index + 1

#######################
# FUNCTIONS
#######################

#test files
INPUT_TEMPLATE = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/input/2INP_core.pdb"
INPUT_QUERY = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/hit_back_bone/PDB_hit_backbone/PDB_search_format5_2INP_core_5tdv-assembly1.cif.gz_F.pdb"

RESNAME_TO_SINGLE_LETTER_DICT = amino_acid_dict = {
    'ALA': 'A',  # Alanine
    'ARG': 'R',  # Arginine
    'ASN': 'N',  # Asparagine
    'ASP': 'D',  # Aspartate
    'CYS': 'C',  # Cysteine
    'GLN': 'Q',  # Glutamine
    'GLU': 'E',  # Glutamate
    'GLY': 'G',  # Glycine
    'HIS': 'H',  # Histidine
    'ILE': 'I',  # Isoleucine
    'LEU': 'L',  # Leucine
    'LYS': 'K',  # Lysine
    'MET': 'M',  # Methionine
    'PHE': 'F',  # Phenylalanine
    'PRO': 'P',  # Proline
    'SER': 'S',  # Serine
    'THR': 'T',  # Threonine
    'TRP': 'W',  # Tryptophan
    'TYR': 'Y',  # Tyrosine
    'VAL': 'V'   # Valine
}

def read_structure_return_sequence(structure: AtomArray) -> str: 
    """
    Read in PDB AtomArray and 
    return single sequence string
    """
    #get all c-alphas from atom array, unique atom per residue
    calpha = structure[structure.atom_name == 'CA']
    
    #iterate over each c_alpha and get string for sequence
    aa_list = [RESNAME_TO_SINGLE_LETTER_DICT[atom.res_name] for atom in calpha]
    aa_seq  = ''.join(letter for letter in aa_list)
    
    return aa_seq

def get_pdb_structure(file: str):
    #get the AtomArray 
    pdb = PDBFile.read(file).get_structure()[0]
    return pdb

def get_ca(structure: AtomArray, type: str = 'start') -> Atom:
    """
    Get either starting or ending CA for a given structure
    """
    if type == 'start':
        start_residue = structure[0].res_id
        print(start_residue, ' is start residue')
        start_ca = structure[(structure.atom_name == 'CA') & (structure.res_id == start_residue)]
        return start_ca

    elif type == 'end':
        end_residue = structure[-1].res_id
        end_ca = structure[(structure.atom_name == 'CA') & (structure.res_id == end_residue)]
        return end_ca
    
    else: 
        print(f'type {type} is not valid')
        exit()

def trim_structure_from_indices(structure: AtomArray, 
                                start_residue: int, 
                                end_residue: int) -> AtomArray:
    """
    Create truncated structure bounded by start and end residue
    """
    assert start_residue < end_residue, f'Start {start_residue} is not less than end {end_residue}'
    return structure[(structure.res_id >= start_residue) & (structure.res_id <= end_residue)]

def find_residue_within_distance(query_atom, hit_pdb):
    """
    Find corresponding start and end position of the template 
    in hit pdb
    """
    shortest = 10
    #create dummy atom
    nearest = Atom([0,0,0])

    for hit_atom in hit_pdb:
        current_distance = distance(hit_atom, query_atom)
        if current_distance < shortest:
            print("current shortest: ", current_distance)
            shortest = current_distance 
            nearest = hit_atom

    return nearest

def write_structure(structure, outFileName):
    """
    Write to output
    """
    outFile = PDBFile()
    outFile.set_structure(structure)
    outFile.write(outFileName)

def trim_ca_based_on_template(input: str, start_ca: Atom, end_ca: Atom) -> AtomArray:
    """
    Find nearest CA and trim structure based on these input atoms
    """
    query_pdb = get_pdb_structure(input)
            
    match_start = find_residue_within_distance(start_ca, query_pdb)
    match_end = find_residue_within_distance(end_ca, query_pdb)

    start_index = match_start.res_id
    end_index = match_end.res_id
        
    print('start and end indices: ', start_index, end_index, )
    trimmed = trim_structure_from_indices(query_pdb, start_index, end_index)

    return trimmed

def superimpose_structure(ref_pdb: AtomArray, alphafold_pdb: AtomArray) -> tuple[AtomArray, float]: 
    """
    Take two PDBFiles and superimpose them
    return the oriented rfdiffusion_pdb and the RMSD
    """
    oriented_alphafold_pdb, transform, fitted_anchor_ind, mobile_anchor_ind = biotite.structure.superimpose_homologs(ref_pdb, alphafold_pdb)
    oriented_alphafold_mask = oriented_alphafold_pdb[mobile_anchor_ind]
    ref_pdb_mask = ref_pdb[fitted_anchor_ind]
    
    rmsd = biotite.structure.rmsd(ref_pdb_mask, oriented_alphafold_mask)
    return oriented_alphafold_pdb, rmsd


def main(): 
    parser = argparse.ArgumentParser("Trimming CA backbone from foldseek output based on template matching")
    parser.add_argument('-t','--template', type=str, help='Template PDB file, the starting structure')
    parser.add_argument('-q','--query', type=str, help='Query PDB file, the hit structure, can be directory of files or single file')
    parser.add_argument('-o', '--outdir', type=str, help='Where to write output file', default='./')

    args = parser.parse_args()

    #do things with the template
    template_pdb = get_pdb_structure(args.template)
    template_basename = os.path.basename(args.template)
    print('template basename', template_basename)
    start_ca = get_ca(template_pdb, 'start')
    end_ca = get_ca(template_pdb, 'end')

    #Check if provided query is single file or directory
    query = args.query
    if os.path.isdir(query): 
        files = [os.path.join(query, file) for file in os.listdir(query) if file.endswith('.pdb')]
        if len(files) == 0:
            print('No .pdb in provided directory')
            exit()

    elif os.path.isfile(query) and query.endswith('.pdb'):
        files = [query]

    else: 
        print(f'{query} is incorrect format')

    #outfilepath
    outfile = os.path.join(args.outdir, 'trim_result.csv')
    os.system(f'echo template,query,trimmed_seq,recalculated_rmsd >> {outfile}')

    for file in files: 
        trimmed = trim_ca_based_on_template(file, start_ca, end_ca)
        
        #get sequence
        seq = read_structure_return_sequence(trimmed)
        
        #do superimposition (likely won't change anything) with trimmed sequence
        aligned_hit, rmsd = superimpose_structure(template_pdb[template_pdb.atom_name == 'CA'], 
                                                  trimmed)
        
        #get basename
        hit_basename = os.path.basename(file).strip('.pdb')
        
        #outpath 
        outpath = os.path.join(args.outdir, f'{hit_basename}_truncated.pdb' )
        write_structure(trimmed, outpath)

        #write to outfile
        os.system(f'echo {template_basename},{hit_basename},{seq},{rmsd} >> {outfile}')
         
if __name__ == '__main__':
    main()
