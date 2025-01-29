# Drug Screening Pipeline Goals

This repository contains scripts for processing Protein Data Bank (PDB) files to analyze amino acid coordinates and their spatial relationship within bounding cuboids. The toolkit provides uses for various geometric operations on proteins, including finding bounding cuboids, placing fixed cuboids, and analyzing amino acid encapsulation and intersection. 

Dependecies:
-biopython
-numpy

## Goals

Goal 1:
For a given pdb file and a set of amino acids determine the 3D coordinates of a bounding cuboid for: 
-Only using the α carbon 
-Using all atoms in the amino acid (aa) 

Goal 2:
For a given pdb, a set of amino acids and a list of lengths (x,y,z lengths of cuboid or a cube if just one length is given) place the fixed cuboid in the center of the list of amino acids for:
-Mean coordinates of α carbon 
-Mean coordinates of all atoms 
-Mean coordinates of side chain atoms 

Goal 3:
For a given pdb and a set of coordinates of a cuboid (or a cube if just one coordinate) find: 
-All the amino acids (and their atoms) fully encapsulated by the cuboid 
-All amino acids α carbons fully encapsulated by the cuboid 
-All amino acids that intersect the cuboid 

For this exercise we will assume that nuclei do not have any volume but they are mere dots.  
Bonus Task: 
For a given region encapsulated by the cuboid find the volume of the space that is not occupied by any atoms (here we will assume that nuclei have volume).  

## Acknowledgments
Special thanks to Alper Celik for reviewing my work and providing guidance.

