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


 Requirements
--------------
 
This code requires you have `make` and `R` on your $PATH and have the
following R packages installed:

  edgeR (tested with 2.6.12)  
  sqldf (tested with 0.4.6.4)  
  gplots (tested with 2.11.0)  
  
R version 2.15.1 was used during development.

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
 

 Details
---------

Currently this is just a placeholder. It will be developed over time
as the paper is finalised.

