BEAM_v1.1.0
Updated November 29, 2021
==================

BEAM was developed by Sudhir Kumar

BEAM is a method aimed at imputing missing bases and correcting base assignment errors in tumor single-cell sequencing data. 

Installation
==================
BEAM is a python script developed in a Windows 64-bit architecture, and does not need to be compiled prior to usage. To use BEAM, the following python packages need to be installed.

Python 3
SciPy
Biopython 
MEGA-CC (please download it from http://www.megasoftware.net/releases/MEGA-CC-7.1.3_win64_setup.exe)

Input file
==================
The input file is the alignment of observed single-cell sequences with MEGA format. Please see Example/Test.meg for an example. 
 
* "T": Mutant allele
* "A": Wild-type allele
* "?": Missing base

Example 
==================
An example dataset (Example/Test.meg) is provided to run BEAM. To run BEAM on Test.meg, please follow commands below from the BEAM directory.

       python BEAM3.py Example/Test.meg

After running BEAM, the output file (refined cell sequences with MEGA format) can be found in the directory of input file. 


How to cite
=================
If you use this BEAM software in your work, please cite the accompanying publication:

Sayaka Miura, Louise A Huuki, Tiffany Buturla, Tracy Vu, and Karen Gomez, and Sudhir Kumar. Computational enhancement of single-cell sequences for inferring tumor evolution. Bioinformatics. 2018 ;34(17):i917-i926 
