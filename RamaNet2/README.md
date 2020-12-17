# RamaNet
Preforms *De novo* protein design using machine learning and PyRosetta to generate a novel synthetic protein structure.

## Requirements:
1. Install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following command (in GNU/Linux) to update your system and install all necessary programs and python modules for this script to run successfully:

`sudo apt update && sudo apt full-upgrade && sudo apt install git dssp gnuplot python3-pip && pip3 install biopython pandas numpy tqdm beautifulsoup4 lxml scipy keras tensorflow`

3. Get this script:

`git clone https://github.com/sarisabban/RamaNet.git`

## Description:
RamaNet (after [Gopalasamudram Ramachandran](https://en.wikipedia.org/wiki/G._N._Ramachandran)) is a script that uses Machine Learning (a GAN network) and PyRosetta to preform *De Novo* Protein Design (from the beginning) i.e. generate and design a synthetic protein structure completely computationally. The script autonomously generates a topology (random every time) then designs a sequence that fits the designed topology using the [RosettaDesign](https://github.com/sarisabban/rosettadesign) protocol. It then submits the structure's FASTA sequence to the [Robetta](http://www.robetta.org/) server to generate and download the custom fragments files in preparation for evaluation using the *Abinitio* folding simulation. The *Abinitio* protocol script can be found [here](https://github.com/sarisabban/RosettaAbinitio). From the fragments files it calculates the RMSD for each fragment's position on the designed structure and plots an (RMSD vs Position) graph to indicate how good the *Abinitio* folding simulation might go (good quality fragments have an average RMSD < 2Å), if the fragments are of bad quality the script will repeat the desgining process up to 10 times before it quits.

## How To Use:
Here is a [video](https://youtu.be/tcHP0IUA7EM) that explains this whole process.
For quick structure generation skip to the last step (step 3).

1. You do not need to generate the Machine Learning dataset, it is already provided and can be downloaded here:

|Dataset name                                                                    |Description                                                                                                                                                                               |
|--------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|[PS+CM dataset](https://www.dropbox.com/s/jmo52kfqy8qf1az/PS%2BCM.tar.bz2?dl=0) |Dataset of all protein structure Φ and Ψ angles and contact map                                                                                                                           |
|[Fragment dataset]()|Dataset of amino acid sequences, secondary structures, SASA, phi, psi, and omega angle features from the vall.jul19.2011 database (used by Rosetta for fragment generation) in .csv format|
|[Sequence dataset]()|Dataset of amino acid sequences !!!! INCOMPLETE |

But if you want to replicate our work use the following command to generate the Machine Learning dataset from the Protein Databank Database, computation time ~168 hours and requires more than 128GB of free disk space.

The default parameters for generating the dataset is isolating proteins between 80 and 150 amino acids, that have more helices and sheets than loops (rigid structures), and with an Rg (radius of gyration) value of less than 15Å (compact structures). The script results in a dataset with the first column as the training example number, then the PDB ID of the file with the chain letter, then the angles *Φ/Ψ/Ω* and the *Cα* distances between carbon-alpha 1 and all other carbon-alphas (hence the PSOC name in the dataset). *0.0* indicates a position with no amino acids, not all protein structures have the same length, but the entire dataset does have the same shape because the empty spaces are filled with zeros. If errors occur, that is fine, some protein files will cause errors (they will be deleted/ignored), but the script should continue all the way to the end and result in a .csv dataset file.

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
* Tabulate all structures' residues' Φ/Ψ/Ω angles and Cα distances.
* Get the largest Cα distance value from the PSOC dataset.

Each of these points can be switched on or off using this strange system we developed. So to run the first section of the protocol use this command:

`python3 RamaNet.py --dataset 111111110000000 DIRECTORY` or `python3 RamaNet.py -d 111111110000000 DIRECTORY`

After the human eye filtering step you can choose what to run, whether to generate the nessesary files to augment the examples through relaxation onto a HPC or just run the augmentation locally, and to generate the Φ/Ψ/Ω angles and Cα distances dataset. The switches are ordered as in the written protocol above.

`python3 RamaNet.py --dataset 000000001001011 DIRECTORY` or `python3 RamaNet.py -d 000000001001011 DIRECTORY`

This system should give you control to run individual step at will.

The most difficult step is the *Human Eye Filtering* step which requires a person to filter out all the unwanted structures manually before moving onto cleaning up each structure and augmenting the data. Unwanted structures such as non-compact structures, structures with more loops than helices and sheets, weird looking structures are all deleted. Also, this is the step to separate structures and collect the ones with traits that you need; the the dataset is augmented (preferably on a HPC to save time). The separation was done manually.

It is best to [contact me](mailto:sari.sabban@gmail.com) if you want to generate your own dataset and I will walk you through the protocol, it is not difficult, but works on individual basis.

2. You do not need to train the neural network because it is already trained and the weights file is available here:

| Weights name                                                                          | Description                                                |
|---------------------------------------------------------------------------------------|------------------------------------------------------------|
|[PSOC Weights]()                                                                       | Neural network weights 

You can use the following command to train the neural network on the dataset (whether you use our dataset or generate your own):

`python3 RamaNet.py --train` or `python3 RamaNet.py -t`

3. Use the following command to generate a novel protein backbone, design the sequence for the backbone, generate fragments from the Robetta server, download them, and analyse the fragment quality:

`python3 RamaNet.py --fragments USERNAME` or `python3 RamaNet.py -f USERNAME`

USERNAME is the username at the Robetta server for fragment generation.

Make sure you have the **weights** file available, either from your training or downloaded from step 2 (provided by us), and that it is in the same directory as the RamaNet2.py script. The file must be named *weights*.

Computation time ~3 hours and will result in 7 files:
* The topology file, which is basically just the structure of the backbone constructed using a sequence of Valines (**backbone.pdb**)
* The final designed structure file - *RosettaDesign* (**structure.pdb**)
* Abinitio input files (**structure.fasta**, **frags.200.3mers**, **frags.200.9mers**, **pre.psipred.ss2**)
* Fragment quality plot (**plot_frag.pdf**)

## References:
When using these scripts kindly reference the following:

* [Sari Sabban, Mikhail Markovsky. (2019) RamaNet: Computational *De Novo* Protein Design using a Long Short-Term Memory Generative Adversarial Neural Network. BioRxiv 671552; doi: https://doi.org/10.1101/671552](https://www.biorxiv.org/content/10.1101/671552v4)

* [Sabban S and Markovsky M. RamaNet: Computational de novo helical protein backbone design using a long short-term memory generative adversarial neural network. F1000Research 2020, 9:298](https://doi.org/10.12688/f1000research.22907.3)
