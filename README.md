 Eczema DRS Paper Code Repository
==================================

Within this repository you will find all the code and data files as 
used/described in the Cole et al paper (see CITATION).

The repository has the following structure:

data/    - source data   
bin/     - scripts for processing the source data  
CITATION - information on how to cite this work  
README   - information on the code and data (this file)  
LICENCE  - information on the licencing of the code  
Makefile - code for running scripts via make  
Makefile.RPTOR - code for running make to do RPTOR-specific analysis

 Requirements
--------------
 
This code requires you have `make` and `R` on your $PATH and have the
following R packages installed:

  edgeR (tested with 3.14.0)  
  sqldf (tested with 0.4.10)  
  gplots (tested with 3.0.1)  
  
R version 3.3.0 was used during development (but should work with anything newer than 2.15.1).

 Quickstart
------------
 
To run all the scripts and generate all figures and data, type:

 `make` 

To run just the analyses and generate data, type:

 `make analysis`

To run generate just the plots, type:

 `make figures`

To delete all outputs, type:

 `make clean`
 
To run RPTOR-specific analysis as per our 2016 paper (see CITATION):

  `make -f Makefile.RPTOR`

 Details
---------

This github repository accompanies two papers and covers the differential gene expression analysis
in a cohort of Irish peadiatric eczema cases in comparison to their filaggrin, FLG, genotype.

The simplest way to use the code is via the Quickstart guide above, but each script found under 
bin/ can run individually if you know what you're doing. Code has some in-line comments.

