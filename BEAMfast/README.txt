BEAMfast_v1.1.0
Updated Feb 8, 2024
==================


BEAM is a method aimed at imputing missing bases and correcting base assignment errors in tumor single-cell sequencing data. 

Installation
==================
BEAM is a python script developed in a Windows 64-bit architecture.
Python 3
SciPy
Biopython 
FastTree.exe (please please place FastTree.exe in this folder)

Input file
==================
The input file is the alignment of observed single-cell sequences with mega format (binary). Please add a germline sequence with the ID of "normal." 

* "T": Mutant allele
* "A": Wild-type allele
* "?": Missing base

How to run
==================

python BEAMfast.py [input mega file]

Output
==================
Output files are binary.

* "T": Mutant allele
* "A": Wild-type allele
* "?": Missing base

How to cite
=================
If you use this BEAM software in your work, please cite the accompanying publication:

Sayaka Miura, Louise A Huuki, Tiffany Buturla, Tracy Vu, and Karen Gomez, and Sudhir Kumar. Computational enhancement of single-cell sequences for inferring tumor evolution. Bioinformatics. 2018 ;34(17):i917-i926 
