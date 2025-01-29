"""
Goal 1:
Module for computing the 3D coordinates of a bounding cuboid for a given PDB file and a set of amino acids.
This can be computed using only the alpha carbons or all atoms in the specified residues.
 
This module uses the BioPython package to parse PDB files and compute the bounding cuboid 
either by only considering alpha carbons or by considering all atoms in specified residues.
Also, numpy is used to perform the necessary calculations.

Dependencies include:
-Biopython
-numpy
"""

#Import packages

from Bio.PDB import *
parser = PDBParser(PERMISSIVE=1) 

import numpy as np

def get_coord_cuboid(pdb_file, amino_acids=None, use_alpha_carbon = True):

    """
    Compute the coordinates of the 3D bounding cuboid for a set of amino acids in a PDB file

    Parameters:
    pdb_file: str
        Path to pdb file to be analyzed.

    amino_acids: list[str], optional, default=None  
        List of amino acid residue names (3-letter codes, e.g, "ALA", "GLY") to be analyzed.
        If None, all residues are analyzed.

    use_alpha_carbon: boolean, default=True
        Whether to use only alpha carbons (Cα) for computing the bounding cuboid.
        If False, all atoms in the residue will be used.

    Output:
    Coordinates of the cuboid:
    -Xmax, Ymax, Zmax: float   
    -Xmin, Ymin, Zmin: float

    Example:
    get_coord_cuboid("example.pdb", amino_acids=["ALA", "GLY"], use_alpha_carbon = True)
    """

    #Load and parse the structure from the PDB file
    structure = parser.get_structure("protein", pdb_file)

    #Initialize an empty list to store the 3D atom coordinates
    coord = []

    #Iterate over the structure (models, chains, residues) to extract the atom coordinates
    for model in structure:
        for chain in model: 
            for residue in chain:
                if amino_acids is None or residue.resname in amino_acids:
                    for atom in residue:
                        if use_alpha_carbon == True:
                            if atom.name == "CA":
                                coord.append(atom.coord)
                        else:
                            coord.append(atom.coord)  

    #Convert the list of coordinates to a numpy array for easier manipulation
    coord_numpy = np.array(coord)

    #Compute bounding cuboid dimensions by finding the maximum and minimum coordinates in each dimension
    x_max,y_max,z_max = np.max(coord_numpy, axis=0) 
    x_min,y_min,z_min = np.min(coord_numpy, axis=0)

    #Print bounding cuboid coordinates 
    print("Xmax=",x_max,"Ymax=",y_max,"Zmax=",z_max,"Xmin=",x_min,"Ymin=",y_min, "Zmin=",z_min)
    