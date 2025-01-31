"""
Goal 3: 
Module for analyzing a given PDB file and a set of coordinates defining a cuboid (or a cube if using just one coordinate) to find: 

a) All the amino acids (and their atoms) fully encapsulated by the cuboid.
b) All amino acids with alpha carbons that are fully encapsulated by the cuboid.
c) All amino acids that intersect the cuboid.
"""

#Import packages
from Bio.PDB import * 
import numpy as np
import pandas as pd

def get_amino_acids_cuboid(pdb_file, cuboid_coord1 = [], cuboid_coord2 = None, cuboid_length = None, use_alpha_carbon = False, intersect = False):

    """
    Extract amino acid and atom information from a PDB file based on spatial constraints defined by a cuboid or cube.

    Parameters:
    pdb_file: str
        Path to pdb file to be analyzed.
    cuboid_coord1: list[float]
        Coordinates [x, y, z] of one corner of the cuboid (or cube center if "cuboid_length" is used).
    cuboid_coord2: list[float], optional, default=None
        Coordinates [x, y, z] of the opposite corner of the cuboid (if not provided, a cube is assumed to be centered at "cuboid_coord1").
    cuboid_length: float, optional, default=None
        Length of the side of the cube 
        Note: if both "cuboid_coord1" and "cuboid_coord2" are provided, "cuboid_length" is ignored.
    use_alpha_carbon: boolean, optional, default=False
        If True, only the alpha carbon (CA) coordinates are used for considered.
    intersect: boolean, optional, default=False
        If True, all amino acids that intersect the cuboid are included.
        Default is false meaning that only fully encapsulated amino acids are considered.

    Output:
    A pandas DataFrame with two columns:
    -"Amino Acid": The residue name of each amino acid.
    -"Atom Type": List of the atom types associated with each amino acid that are encapsulated or intersect with the cuboid.

    Notes
    - If "use_alpha_carbon" is set to True, only alpha carbon atoms are checked against the cuboid bounds.
    - When "intersect=True", any amino acid with at least one atom inside the cuboid is included.
    - Fully encapsulated amino acids require all atoms to fall within the cuboid boundaries.

    Example: 
    get_amino_acids_cuboid("example.pdb", cuboid_coord1 = [12,14,8], cuboid_coord2 = [-12,-12,-12], use_alpha_carbon = False, intersect = True)
    """

    #Check that the coordinates are valid
    if len(cuboid_coord1) != 3:
        raise ValueError("Invalid coordinates entered: invalid number of coordinates. Must enter x,y,z coordinates for cuboid_coord1")
    
    #Check if the cuboid is a cuboid or a cube and set the max and min coordinates accordingly
    if cuboid_coord2 is not None:
        if all(x > y for x, y in zip(cuboid_coord1, cuboid_coord2)):
            x_max,y_max,z_max = cuboid_coord1
            x_min,y_min,z_min = cuboid_coord2
        elif all(x < y for x, y in zip(cuboid_coord1, cuboid_coord2)):
            x_max,y_max,z_max = cuboid_coord2
            x_min,y_min,z_min = cuboid_coord1
        else:
            raise ValueError("Invalid coordinates entered: inconsistent dimensions")
    elif cuboid_coord2 is None and cuboid_length is not None:
        x_max = cuboid_coord1[0] + cuboid_length/2
        x_min = cuboid_coord1[0] - cuboid_length/2
        y_max = cuboid_coord1[1] + cuboid_length/2  
        y_min = cuboid_coord1[1] - cuboid_length/2
        z_max = cuboid_coord1[2] + cuboid_length/2
        z_min = cuboid_coord1[2] - cuboid_length/2
    elif cuboid_coord2 is None and cuboid_length is None:
        raise ValueError("Invalid coordinates entered: must provide either cuboid_coord2 or cuboid_length")

    #Load and parse the structure from the PDB file
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("protein", pdb_file)

    #Initialize empty lists to store the 3D atom coordinates and the atom type
    amino_acid_list = []
    atom_type_list = []

    #Iterate over the structure (models, chains, residues) to extract the atom coordinates
    #Check if the atoms are within the cuboid boundaries and store the amino acid and atom type information according to the specified conditions
    for model in structure:
        for chain in model: 
            for residue in chain:
                residue_coord = []
                atom_type = []
                for atom in residue:
                    if use_alpha_carbon == True:
                        if atom.name == "CA":
                            atom_coord = atom.coord
                            if (x_max >= atom_coord[0] >= x_min and y_max >= atom_coord[1] >= y_min and z_max >= atom_coord[2] >= z_min):
                                amino_acid_list.append(residue.resname)
                                atom_type_list.append(atom.name)
                    elif use_alpha_carbon == False and intersect == True:
                        residue_coord = [atom.coord for atom in residue]
                        atom_type = [atom.name for atom in residue]
                        if any(x_max >= atom_coord[0] >= x_min and y_max >= atom_coord[1] >= y_min and z_max >= atom_coord[2] >= z_min for atom_coord in residue_coord):
                            amino_acid_list.append(residue.resname)
                            atom_type_list.append(tuple(atom_type))
                    elif use_alpha_carbon == False and intersect == False:
                        residue_coord = [atom.coord for atom in residue]
                        atom_type = [atom.name for atom in residue]
                        if all(x_max >= atom_coord[0] >= x_min and y_max >= atom_coord[1] >= y_min and z_max >= atom_coord[2] >= z_min for atom_coord in residue_coord):
                            amino_acid_list.append(residue.resname)
                            atom_type_list.extend(atom_type)

    #Create a pandas dataframe with the amino acid and atom type information                           
    data = (list(zip(amino_acid_list, atom_type_list)))
    df = pd.DataFrame(data, columns = ["Amino Acid", "Atom Type"])

    #Remove any duplicate amino acids and atoms
    df['Amino Acid'] = df['Amino Acid'].str.strip()
    df = df.drop_duplicates()

    #Remove any water or zinc atoms (non amino acid atoms)
    df = df[~df["Amino Acid"].isin(["ZN", "HOH"])]

    #Group the atoms by amino acid
    df = df.groupby("Amino Acid")["Atom Type"].apply(list).reset_index()

    #Return dataframe
    return(df)
    
                       




