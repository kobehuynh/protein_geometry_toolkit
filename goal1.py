#Goal 1:For a given pdb file and a set of amino acids determine the 3D coordinates of a bounding cuboid for  
#a) Only using the α carbon 
#b) Using all atoms in the amino acid (aa) 

#Import packages

from Bio.PDB import *
parser = PDBParser(PERMISSIVE=1) 

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#Create function to obtain x,y,z coordinates of bounding cuboid given a pdb file and creates a visualization of the 3D bounding cuboid
#Can set list of amino acids or determine whether alpha carbons should be included

def get_coord_cuboid(pdb_file, amino_acids=None, use_alpha_carbon = True):
    
    structure = parser.get_structure("protein", pdb_file)

    coord = []

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

    coord_numpy = np.array(coord)

    x_max,y_max,z_max = np.max(coord_numpy, axis=0) 
    x_min,y_min,z_min = np.min(coord_numpy, axis=0)

    vertices = [
        [x_min, y_min, z_min],  
        [x_min, y_max, z_min],
        [x_max, y_max, z_min],
        [x_max, y_min, z_min],
        [x_min, y_min, z_max],  
        [x_min, y_max, z_max],
        [x_max, y_max, z_max],
        [x_max, y_min, z_max],
    ]
    faces = [
        [vertices[0], vertices[1], vertices[2], vertices[3]],  
        [vertices[4], vertices[5], vertices[6], vertices[7]],  
        [vertices[0], vertices[1], vertices[5], vertices[4]],  
        [vertices[2], vertices[3], vertices[7], vertices[6]], 
        [vertices[1], vertices[2], vertices[6], vertices[5]],  
        [vertices[0], vertices[3], vertices[7], vertices[4]],  
    ]

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    poly3d = Poly3DCollection(faces, alpha=0.3, facecolors = "green", edgecolors = "black")
    ax.add_collection3d(poly3d)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Visualization of 3D Bounding Cuboid")
    plt.show()

    print("Xmax=",x_max,"Ymax=",y_max,"Zmax=",z_max,"Xmin=",x_min,"Ymin=",y_min, "Zmin=",z_min)
    
