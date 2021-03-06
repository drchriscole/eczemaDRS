
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3903510.svg)](https://doi.org/10.5281/zenodo.3903510)


Eczema DRS Paper Code Repository
==================================

Within this repository you will find all the code and data files as 
used/described in the Cole et al paper (see CITATION).

It has also been adapted for work on a collaboration with Ryan O'Shaughnessy 
(UCL) on RPTOR gene expression (see CITATION).

The repository has the following structure:

app/      - Shiny app for visualising gene expression
bin/      - scripts for processing the source data  
data/     - source data   
CITATION  - information on how to cite this work  
README.md - information on the code and data (this file)  
LICENCE   - information on the licencing of the code  
Makefile  - code for running scripts via make  
Makefile.RPTOR - code for performing RPTOR-specific analysis

 Requirements
--------------
 
This code requires you have `make` and `R` on your $PATH and have the
following R packages installed:

  edgeR (tested with 3.14.0)  
  sqldf (tested with 0.4.10)  
  gplots (tested with 3.0.1)  

For the Shiny App the `shiny` R library will also be required.
  
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
bin/ can run be individually if you know what you're doing. Code has some in-line comments.

 Shiny App
-----------

A utility webapp has also been developed with RStudio's Shiny framework to perform gene-specific searches of expression profiles between cases and controls or stratified by FLG genotype.

To run the Shiny App do the following:

  1. Clone this repository into [a new RStudio project](https://happygitwithr.com/rstudio-git-github.html#clone-the-new-github-repository-to-your-computer-via-rstudio)
  2. Open the `app\DRSexpr\app.R` file
  3. Click the 'Run App' button

