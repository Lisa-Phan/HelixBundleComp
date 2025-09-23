"""
Lisa P. 
9/17/25

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
"""
import argparse
import biotite
import biotite.sequence as seq
import biotite.sequence.align as align
from biotite.structure.io.pdb import PDBFile
from biotite.structure import distance
from biotite.structure import Atom


#test files
INPUT_QUERY = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/input/2INP_core.pdb"
INPUT_TEST = r"/work/09069/dhp563/ls6/HelixBundleComp/structure_query/hit_back_bone/PDB_hit_backbone/PDB_search_format5_2INP_core_5tdv-assembly1.cif.gz_F.pdb"

def get_pdb_structure(file: str):
    #get the AtomArray 
    pdb = PDBFile.read(file).get_structure()[0]
    return pdb

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

def read_pdb_return_sequence(structure):
    pass



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
    
    print(alignment)





def test():
# test function to parse 

    query_pdb = get_pdb_structure(INPUT_QUERY)
    hit_pdb = get_pdb_structure(INPUT_TEST)
        
    start_atom = query_pdb[0]
    end_atom = query_pdb[-1]

    print(vars(start_atom))

    # match_start = find_residue_within_distance(start_atom, hit_pdb)
    # match_end = find_residue_within_distance(end_atom, hit_pdb)

    # print(match_start)
    # print(match_end)


if __name__ == '__main__':
    test()
