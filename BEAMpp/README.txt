BEAMpp_v1.1.0
Updated Feb 3, 2024
==================

The posterior probability of observing the mutant base at each position of a cell will be calculated for a given cell sequence alighment and tree. 

Installation
==================
BEAM is a python script developed in a Windows 64-bit architecture. To use BEAM, the following python packages need to be installed.

Python 3
SciPy
Biopython 

Input file
==================
1. Cell sequence alignment.
The input file is the alignment of observed single-cell sequences with fasta format, and the extension should be .fasta. A germline cell sequence needs to be included. Please see Example/Test.fasta for an example. 

2. Cell phylogeny
Please infer a cell phylogeny using the input fasta file and place it in the same folder as the fasta file. The file ID needs to be the same as the fasta file, and the extension should be .nwk. e.g., Example/Test.nwk. 

How to run
==================

python BEAMpp.py [input fasta file] [germline cell ID]

Example 
==================

python BEAMpp.py Example/Test.fasta normal

After running BEAM, the output file (refined cell sequences) can be found in the directory of input file. 

A: Non-mutant base
T: Mutant base

The posterior probability is listed in the text file. 

How to cite
=================
If you use this BEAM software in your work, please cite the accompanying publication:

Sayaka Miura, Louise A Huuki, Tiffany Buturla, Tracy Vu, and Karen Gomez, and Sudhir Kumar. Computational enhancement of single-cell sequences for inferring tumor evolution. Bioinformatics. 2018 ;34(17):i917-i926 
