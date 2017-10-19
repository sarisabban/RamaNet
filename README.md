# ResottaDeNovo
Preforms De Novo Design using Machine Learning and PyRosetta to generates a novel protein structure.

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following command (in GNU/Linux) to install all necessary programs and python modules for this script to run successfully:

`sudo apt update && sudo apt install python3-pip DSSP gnuplot && sudo pip3 install biopython bs4`

## Description:
This is a script that uses Machine Learning and PyRosetta to preform De Novo Design (from the beginning) i.e. develop and design a synthetic protein structure totally computationally. There is no input for this script, it autonomously generates a topology (random every time) then designs a sequence that fits this topology, then submits the structure's FASTA sequence to the [Robetta](http://www.robetta.org/) server to generate and download the custom fragment files in preparation for an Abinitio folding simulation. Abinitio folding script can be found [here](https://github.com/sarisabban/RosettaAbinitio). Finally it calculates the RMSD for each fragment position on the designed structure and plots an (RMSD vs Position) graph to indicate how good the Abinitio folding simulation might go (idealy you want all positions to be under 2Å RMSD).

The script will generate 1 structure1. It is advised to run this script in an array to generate multiple strctures and see which one has low RMSD fragments. Mind you, if you generate too many structures this might overwhelm the Robetta Server by submitting and requesting too many fragment files, please be considirate and run this script once to generate 1 structure at a time only.

## How To Use:
1. Use the following command to run the script:

`python3 DeNovo.py`

2. Computation time is very long (around 12 hours).
3. The Script will result in 7 files:
* Topology file (DeNovo.pdb)
* Sequence designed file (structure.pdb)
* Abinitio input files (structure.pdb, frags.200.3mers, frags.200.9mers, pre.psipred.ss2)
* Fragment quality plot (plot_frag.pdb)

# This script is still under development. This statement will be removed when the script is complete and benchmarked.
