# ResottaDeNovo
Preforms De Novo Design using Machine Learning and PyRosetta to generates a novel synthesic protein structure.

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following command (in GNU/Linux) to install all necessary programs and python modules for this script to run successfully:

`sudo apt update && sudo apt install python3-pip DSSP gnuplot && sudo pip3 install biopython biopandas bs4 tqdm`

## Description:
This is a script that uses Machine Learning and PyRosetta to preform De Novo Design (from the beginning) i.e. develop and design a synthetic protein structure totally computationally. There is no input for this script, it autonomously generates a topology (random every time) then designs a sequence that fits this topology, then submits the structure's FASTA sequence to the [Robetta](http://www.robetta.org/) server to generate and download the custom fragment files in preparation for an Abinitio fold simulation. The Abinitio script can be found [here](https://github.com/sarisabban/RosettaAbinitio). Finally it calculates the RMSD for each fragment position on the designed structure and plots an (RMSD vs Position) graph to indicate how good the Abinitio fold simulation might go (idealy you want all positions to be under 2Å RMSD).

The script will generate 1 structure. It is advised to run this script in an array to generate multiple strctures and see which one has low RMSD fragments. Mind you, if you generate too many structures this might overwhelm the Robetta Server by submitting and requesting too many fragment files, please be considirate and run this script once to generate 1 structure at a time only.

## How To Use:
1. Use the following command to generate the Machine Learning Dataset from the Protein Databank Database (computation time ~72 hours and requires ~120GB of free disk space):

`python3 Database.py`

This script will result in 1 file:

* data.csv

The default parameters for the Database.py script is isolating proteins between 100 and 150 amino acids, and with an Rg value of less than 15. The script results in a dataset with the first column is the training example number, then the name of the PDB file, then the type of secondary structure at each amino acid position [0 = Nothing or empty (i.e a short protein) , 1 = Loop , 2 = Helix , 3 = Strand], then there are the constraint distances starting with the distance in angstrom between the first and last amino acids, then between 10 amino acids from the first position and 10 amino acids from the last position, etc...

2. Use the following command to preform DeNovo protein design:

`python3 DeNovo.py`

This script will result in 7 files:
* Topology file, which is basically just the structure of the backbone drawn using a sequence of Valine (DeNovo.pdb)
* Sequence designed file (structure.pdb)
* Abinitio input files (structure.pdb, frags.200.3mers, frags.200.9mers, pre.psipred.ss2)
* Fragment quality plot (plot_frag.pdb)

# This script is still under development. This statement will be removed when the script is completed and benchmarked.
