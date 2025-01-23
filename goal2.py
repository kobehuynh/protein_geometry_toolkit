#Goal 2: For a given pdb, a set of amino acids and a list of lengths (x,y,z lengths of cuboid or a cube 
#if just one length is given) place the fixed cuboid in the center of the list of amino acids for 

#a) Mean coordinates of α carbon 
#b) Mean coordinates of all atoms 
#c) Mean coordinates of side chain atoms 

#Import packages

from Bio.PDB import * 
import numpy as np

#Function to calculate the coordinates of a cuboid given a pdb file, a list of amino acids, a mean coordinate type 
#(alpha_carbons, all_atoms, side_chain_atoms) and a list of cuboid dimensions (x,y,z lengths of cuboid or a cube if just one length is given)

def get_coord_cuboid_center(pdb_file, amino_acids=None, mean_coord_type = "all_atoms", cuboid_dimensions = []):
    
    valid_coord_types = ["alpha_carbons","all_atoms","side_chain_atoms"]
    if mean_coord_type not in valid_coord_types:
        raise ValueError("Invalid coordinate type entered. Please enter one of the following:" "alpha_carbons", "all_atoms", "side_chain_atoms")

    if cuboid_dimensions is None:
        raise ValueError("You must enter values for cuboid dimensions")

    #If only one value is entered, assume that the cuboid is a cube
    if len(cuboid_dimensions) == 1:
        cuboid_dimensions = np.array([cuboid_dimensions[0]] * 3)
    elif len(cuboid_dimensions) != 3:
        raise ValueError("Incorrect number of values for cuboid dimensions")
    else:
        cuboid_dimensions = np.array(cuboid_dimensions)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("protein", pdb_file)

    coord = []

    for model in structure:
        for chain in model: 
            for residue in chain:
                if amino_acids is None or residue.resname in amino_acids:
                    for atom in residue:
                        if mean_coord_type == "alpha_carbons":
                            if atom.name == "CA":
                                coord.append(atom.coord)
                        elif mean_coord_type == "all_atoms":
                            coord.append(atom.coord)  
                        elif mean_coord_type == "side_chain_atoms": 
                            if atom.name not in ("CA","O","N","C"):
                                coord.append(atom.coord)

    coord_numpy = np.array(coord)

    x_mean,y_mean,z_mean = np.mean(coord_numpy, axis=0) 

    x_max = x_mean + cuboid_dimensions[0]/2
    x_min = x_mean - cuboid_dimensions[0]/2
    y_max = y_mean + cuboid_dimensions[1]/2
    y_min = y_mean - cuboid_dimensions[1]/2
    z_max = z_mean + cuboid_dimensions[2]/2
    z_min = z_mean - cuboid_dimensions[2]/2

    return("Xmax=",x_max,"Ymax=",y_max,"Zmax=",z_max,"Xmin=",x_min,"Ymin=",y_min, "Zmin=",z_min)
    
