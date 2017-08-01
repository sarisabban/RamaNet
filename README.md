# ResottaDeNovo
Preforms De Novo Design using PyRosetta. Generates a novel protein with only helices (no sheet).

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following commands (in GNU/Linux) to install all nessesary Python libraries for this script to run successfully:

`sudo apt install python3-pip pymol DSSP gnuplot`

`sudo python3 -m pip install numpy biopython bs4`

## Description:
This is a script that uses PyRosetta to preform De Novo Design (from the beginning) i.e. develop and design a synthetic protein structure totally computationally. There is no input for this script, it autonomously generates a topology (random every time) then designs a sequence that fits this topology, then submits the structure's FASTA sequence to the [Robetta](http://www.robetta.org/) server to generate and download the custom fragment files in preparation for an Abinitio folding simulation. Abinitio folding script can be found [here](https://github.com/sarisabban/RosettaAbinitio). Finally it calculates the RMSD for each fragment position on the designed structure and plots a (RMSD vs Position) graph to indicate good the Abinitio folding simulation might go (idealy you want all positions to be under 1Å RMSD).

The script will generate 1 structure. It is advised to run this script in an array to generate multiple strctures and see which one has low RMSD fragments. Mind you, if you generate too much structures this might overwhelm the Robetta Server but submitting and requesting too many fragment files, please be considirate and generate 2 structures at a time only.

## How To Use:
1. Use the following command to run the script:

`python3 DeNovo.py`

2. Calculation time is very long.
3. The Script will result in 7 files:
* Topology file (DeNovo.pdb)
* Sequence designed file (structure.pdb)
* Abinitio input files (structure.pdb, t000_.psipred_ss2, aat000_03_05.200_v1_3, aat000_09_05.200_v1_3)
* Fragment quality plot (plot_frag.pdb)

# This script is still under development. If you can contribute please do, and your name will be included.
