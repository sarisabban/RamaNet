# RamaNet
Preforms *De novo* protein design using machine learning and PyRosetta to generate a novel synthetic protein structure.

## Requirements:
1. Install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following command (in GNU/Linux) to update your system and install all necessary programs and python modules for this script to run successfully:

`sudo apt update && sudo apt full-upgrade && sudo apt install dssp gnuplot python3-pip && pip3 install biopython pandas numpy matplotlib tqdm beautifulsoup4 lxml scipy keras tensorflow`

## Description:
RamaNet (after [Gopalasamudram Ramachandran](https://en.wikipedia.org/wiki/G._N._Ramachandran)) is a script that uses Machine Learning (a GAN network) and PyRosetta to preform *De Novo* Protein Design (from the beginning) i.e. generate and design a synthetic protein structure completely computationally. There is no input for this script, it autonomously generates a topology (random every time) then designs a sequence that fits the designed topology using the [RosettaDesign](https://github.com/sarisabban/rosettadesign) protocol. It then submits the structure's FASTA sequence to the [Robetta](http://www.robetta.org/) server to generate and download the custom fragments files in preparation for an *Abinitio* folding simulation. The *Abinitio* protocol script can be found [here](https://github.com/sarisabban/RosettaAbinitio). Finally it calculates the RMSD for each fragment's position on the designed structure and plots an (RMSD vs Position) graph to indicate how good the *Abinitio* folding simulation might go (ideally you want an average RMSD < 2Å).

The script will generate one structure. It is advised to run this script and generate multiple structures and see which one has the lowest average RMSD fragments. But mind you, if you generate too many structures this might overwhelm the Robetta server by submitting and requesting too many fragment files, please be considerate and run this script once to generate one structure at a time only.

## How To Use:
For a quick structure generation right now just skip to the last step (step 4).

1. You do not need to generate the Machine Learning datasets, they are already provided and can be downloaded here:

| Dataset name                                                                        | Description                                                                |
|-------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
|[Helix PS dataset](https://www.dropbox.com/s/wdi7dxmshgwexuk/PS_Helix_500.csv?dl=0)  | Dataset of only helical structures' Phi/Psi angles                         |
|[Helix PSC dataset](https://www.dropbox.com/s/3mg6edh933uhzu8/PSC_Helix_500.csv?dl=0)| Dataset of only helical structures' Phi/Psi angles and Constraints         |
|[Sheet PSC dataset](https://www.dropbox.com/s/ws1zelxl2jm1n3j/PSC_Sheet_500.csv?dl=0)| Dataset of only sheet structures' Phi/Psi angles and Constraints           |
|[Mix PSC dataset](https://www.dropbox.com/s/qz35dsgvs91wsjz/PSC_Mix_500.csv?dl=0)    | Dataset of mixed helix and sheet structures' Phi/Psi angles and Constraints|

But if you want to replicate our work use the following command to generate the Machine Learning dataset from the Protein Databank Database (computation time ~168 hours and requires more than 128GB of free disk space):

`python3 Database.py`

The default parameters for the Database.py script is isolating proteins between 80 and 150 amino acids, that have more helices and strands than loops (a rigid structure), and with an Rg value of less than 15 (compact structure). The script results in a dataset with the first column as the training example number, then the PDB ID of the file (and chain letter), then the angles *Phi/Psi* for each amino acid (PS dataset), an option is to also include the cst constraint values between carbon-alpha 1 and all other carbon-alphas (PSC dataset). *0.0* indicates a position with no amino acids, not all protein structures have the same length, but the entire dataset does have the same length and shape because the empty spaces are filled with zeros. If errors occur, that is fine, some protein files will cause errors (and they will be deleted/ignored), but the script should continue all the way to the end and result in a dataset file. 

The dataset generation protocol is as follows:
* Download the PDB database
* Extract files
* Remove non-protein structures
* Remove structures less than or larger than a specified amino acid length
* Remove structure with broken chains
* Remove structures that have loops that are larger than a specific length
* Renumber structures' amino acids
* Remove structures that are below a specified Radius of Gyration value
* ########## --- HUMAN EYE FILTERING --- ##########
* Clean every structure in the database
* Make a list of all paths (if next step is performed in a high performance computer HPC)
* Generate HPC submission file (PBS job scheduler)
* Relax each structure multiple times (on a HPC or a local computer), this is to augment the examples
* Get each residue's phi/psi angles (and constraints if required - you have to comment in and out the relevent functions at the end of the script)

The most difficult step is the *Human Eye Filtering* step which requires a person to filter out all the unwanted structures manually before moving onto cleaning up each structure and augmenting the data which eventually results in a .csv file. Unwanted structures such as non-compact structures, structures with more loops than helices and sheets, weird looking structures. Also, this is the step to separate structures and collect the ones with traits that you need; I decided to separate the dataset into structures with only helices, only sheet, and a mix of the two before augmenting each dataset. This was in order to walk the neural network slowly through the training process. The separation was done manually.

It is best to [contact me](mailto:sari.sabban@gmail.com) if you want to generate your own database and I will walk you through the protocol, it is not difficult, but works on individual basis. If you contact me I can also identify what people do not understand about this particular script and I will be able to modify this README file to make it easier to use.

2. You do not need to train the neural neural network because it is already trained and the weights file is available here:

| Weights name                                                                          | Description                                                |
|---------------------------------------------------------------------------------------|------------------------------------------------------------|
|[Helix PS Weights](https://www.dropbox.com/s/ojv1ugryj4tqpnm/PS_Helix_Weights.zip?dl=0)| Neural network weights generated from the Helix PS dataset |
|[Helix PSC Weights]()                                                                  | Neural network weights generated from the Helix PSC dataset|
|[Sheet PSC Weights]()                                                                  | Neural network weights generated from the Sheet PSC dataset|
|[Mix PSC Weights]()                                                                    | Neural network weights generated from the Mix PSC dataset  |

You can use the following command to train the neural network on the dataset (whether you use ours or generate your own):

`python3 Generate.py --train` or `python3 Generate.py -t`

3. Use the following command to generate a novel protein structure, generate fragments from the Robetta server, download these fragment files, and analyse these fragments:

`python3 Generate.py --fragments USERNAME`

USERNAME is the username at the Robetta server for fragment generation.

4. Use the following command to only generate a novel protein structure without generating any fragments:

`python3 Generate.py`

Make sure you have the **weights** directory available, either from your training or downloaded from step 2 (provided by us), and that it is in the same directory as the Generate.py script. The directory must be named *weights*.

This script (computation time ~24 hours) will result in 1 file (if no fragments are requested), or 7 files (if fragments are requested):
* Topology file, which is basically just the structure of the backbone drawn using a sequence of Valine (**backbone.pdb**)
* The final designed structure file - *RosettaDesign* (**structure.pdb**)
* Abinitio input files (**structure.fasta**, **frags.200.3mers**, **frags.200.9mers**, **pre.psipred.ss2**)
* Fragment quality plot (**plot_frag.pdf**)

## References:
When using these scripts kindly reference the following:
* 
