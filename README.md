# PASM- phylogenetic analysis-MSA statistics-mutation identification

**PASM is a python program that builds a phylogenetic tree, calculates multiple sequence alignment statistics, and identifies key mutations. It performs Needleman–Wunsch global alignment, builds the distance matrix with the kimura-2 [k80] model and constructs the phylogenetic tree using the UPGMA algorithm. The conservation score and the key mutations identified from the multiple sequence alignment file are calculated, and the accuracy of the phylogenetic tree is validated by implementing the bootstrap method. The program is useful for biologists and bioinformaticians in viral evolution studies.**

# Required Python packages
**Biopython, numpy, matplotlib**

# Algorithms :
Global alignment: Needleman–Wunsch
Distance model: Kimura 2-Parameter (K80/K2P)
Tree construction: Unweighted Pair Group Method with Arithmetic Mean (UPGMA) algorithm
Conservation scoring: Shannon entropy
Tree validation: Bootstrap analysis

# How to run the program:
**Input: The user can input a FASTA file of sequences and use the file_to_folder() function to convert the file to a folder of FASTA files of one sequence or directly use a folder as input.**
In Line 13
path = "specify path to save the folder"
In Line 42
file_path = "enter the fasta file path"
In Line 43
file_to_folder(file_path,name of the folder)In line 
In line 130:
os.chdir("enter the path of where you want to save the pairwise alignment output files")
In line 201: 
dataset_path ="enter the path of the folder that contain the sequences"
In line 647:
for i in os.listdir("enter the path of the folder that contain the sequences"):
In line 656:
stats=Statistics("enter the path of the multiple sequence alignment file to calculate the alignment statistics (NAdata_MSAoutput2021.FASTA) or ()")

# Project structure
**The program is composed of 3 main classes 
Global alignment class (calculate the distance matrix)
UPGMA class (construct phylogenetic tree)
Statistics class (calculate all statistics, identify mutation, validate the tree [bootstrap], visualize the tree)**

# Output
**phylogenetic tree, text file for alignment statistics, text file for the key mutation, FASTA file for the pairwise alignment if needed**
