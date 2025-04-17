"""
Goal 1:
This module computes the 3D bounding cuboid coordinates that encapsulate a specified set of residues
from a given PDB file. The cuboid is computed using either all atoms or only the alpha carbon (Cα) atoms
of the specified residues.
 
This module uses the BioPython package to parse PDB files and compute the bounding cuboid 
either by only considering alpha carbons or by considering all atoms in specified residues.
Also, numpy is used to perform the necessary calculations.

Dependencies include:
-Biopython
-numpy
"""

from Bio.PDB import PDBParser
import numpy as np

def get_coord_cuboid(pdb_file, chain_id, residue_number, use_alpha_carbon = True):

    """
    Compute the coordinates of the 3D bounding cuboid for a set of amino acids in a PDB file

    Parameters:
    pdb_file: str
        Path to pdb file to be analyzed.

    chain_id : list of str
        List of chain identifiers corresponding to the target residues.

    residue_number : list of int
        List of residue numbers corresponding to the target residues.
        Must be the same length and order as `chain_id`.

    use_alpha_carbon: boolean, default=True
        Whether to use only alpha carbons (Cα) for computing the bounding cuboid.
        If False, all atoms in the residue will be used.
    
    Output:
    Coordinates of the cuboid:
    -Xmax, Ymax, Zmax: float   
    -Xmin, Ymin, Zmin: float

    Example:
    get_coord_cuboid("your_file.pdb", chain_id = ["A","A","A","A","A","B","B","C"], residue_number= [1,2,3,4,5,4,5,7], use_alpha_carbon = True)
    """

    if len(chain_id) != len(residue_number):
        raise ValueError("'chain_id' list and 'residue_number' list must be of the same length")

    #Create list of tuples using chain ids and residue numbers
    residue_tuple = list(zip(chain_id, residue_number))

    #Load and parse the structure from the PDB file
    parser = PDBParser(PERMISSIVE=1) 
    structure = parser.get_structure("protein", pdb_file)

    #Initialize an empty list to store the 3D atom coordinates
    coord = []

    #Iterate over the structure (models, chains, residues) to extract the atom coordinates
    for model in structure:
        for chain in model: 
            for residue in chain:
                hetfield, resseq, icode = residue.get_id()
                if hetfield == ' ' and (chain.id, resseq) in residue_tuple:
                    for atom in residue:         
                        if not use_alpha_carbon or atom.name == "CA":
                            coord.append(atom.coord)
            
    if not coord:
        raise ValueError("No coordinates were collected. Check your chain/residue input or PDB file.")

    #Convert the list of coordinates to a numpy array for easier manipulation
    coord_numpy = np.array(coord)

    #Compute bounding cuboid dimensions by finding the maximum and minimum coordinates in each dimension
    x_max,y_max,z_max = np.max(coord_numpy, axis=0) 
    x_min,y_min,z_min = np.min(coord_numpy, axis=0)

    return {
        "Xmax": x_max,
        "Ymax": y_max,
        "Zmax": z_max,
        "Xmin": x_min,
        "Ymin": y_min,
        "Zmin": z_min
    }


