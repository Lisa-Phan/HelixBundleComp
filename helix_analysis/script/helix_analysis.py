"""
8/9/25
Lisa P. 

Test mda implementation of HELANAL
HELANAL needs residue specifications, so first need to find 
a way to pick out which residue is part of the helix with different algorithm first

TODO: write a generic function that takes in PDB file
run DSSP on it to assign helix, and then return residue selection
that specifies start and end of helix resnum pairs
"""

import MDAnalysis as mda
from MDAnalysis.analysis import helix_analysis as hel
from MDAnalysis.analysis.dssp import DSSP

RNR_PDB=r"/stor/scratch/YiLu/dhp563/HelixBundleComp/structure_query/input/1W69_core.pdb"
PHP_PDB=r"/stor/scratch/YiLu/dhp563/HelixBundleComp/structure_query/input/2INP_core.pdb"
TMO_PDB=r"/stor/scratch/YiLu/dhp563/HelixBundleComp/structure_query/input/5TDU_core.pdb"
MMO_PDB=r"/stor/scratch/YiLu/dhp563/HelixBundleComp/structure_query/input/6YD0_core.pdb"

def record_local_bend_angles(structure: str, resnum_tuple: tuple) -> dict:
    """
    structure: string for path of PDB file
    resnum_tuple, a tuple of tuples, each of which is two integers specifying 
    the start and end of a helix

    return a dictionary of key: startindex_endindex : [local bend angles]
    """ 
    local_bend_angle_dict = {}
    universe = mda.Universe(structure)
    for value in resnum_tuple:
        start_resn, end_resn = value
        key = f'{start_resn}-{end_resn}'
        helanal = hel.HELANAL(universe, select=f'name CA and resnum {key}')
        helanal.run()
        print(helanal.results.local_bends)
        local_bend_angle_dict[key] = helanal.results.local_bends[0]
    return local_bend_angle_dict

RNR = record_local_bend_angles(RNR_PDB, ((121, 151), (222, 250), (253, 280)))
PHP = record_local_bend_angles(PHP_PDB, ((91, 120), (124, 156), (170, 184), (188, 216), (220, 249)))
TMO = record_local_bend_angles(TMO_PDB, ((87, 117), (121, 163), (167, 181), (185, 214), (218, 247)))
MMO = record_local_bend_angles(MMO_PDB, ((97, 126), (131, 161), (174, 193), (197, 227), (231, 257)))

print(RNR, PHP, TMO, MMO)
