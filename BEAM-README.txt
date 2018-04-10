CloneScape_v0.1.0
Updated April 10, 2018
==================

BEAM was developed by Sudhir Kumar

BEAM is a method aimed at imputing missing bases and correcting base assignment errors in tumor single-cell sequencing data. 

Installation
==================
BEAM is a python script developed in a Windows and Unix/Linux 64-bit architecture, and does not need to be compiled prior to usage. To use BEAM, the following python packages need to be installed.

Python 2
SciPy
Biopython 

Additionally, BEAM uses MEGA-CC (>= 7.0.18). MEGA-CC can be downloaded for Windows, Mac OS X, and Linux from (http://www.megasoftware.net).

Input file
==================
The input file is the alignment of observed single-cell sequences with MEGA format. Please see Example/Test.meg for an example. 
 
* "T": Mutant allele
* "A": Wild-type allele
* "?": Missing base

Example 
==================
An example dataset (Example/Test.meg) is provided to run BEAM. To run BEAM on Test.meg, please follow commands below from the BEAM directory.

       python BEAM.py path/Test.meg

After running BEAM, the output file (refined cell sequences with MEGA format) can be found in the directory of input file. 


How to cite
=================
If you use this BEAM software in your work, please cite the accompanying publication:

Sayaka Miura, Louise A Huuki, Tiffany Buturla, Tracy Vu, and Karen Gomez, and Sudhir Kumar. Imputing missing data and correcting base assignment errors in tumor single cell sequences. XXX (2018) 
