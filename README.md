# Scripts
Scripts used for analysis of genomic data. 

This repos contains modular scripts for the following dataset types:
* ChIP-seq
* CUT&RUN
* RNA-seq
* ATAC-seq
* WGS & polymorphic mobile element insertion calling

These scripts are designed to work on CU Boulder's computing system, Fiji, which uses Slurm as a resource manager and job scheduler. 

These scripts are intentionally modular - each should contain constant variables that, ideally, should not need to be changed between projects. Project specific path variables, such as input and output directories, will be specified at the command line when submitting the script.

Finally, these scripts are a work in progress. Some time in the future I intend to overhaul these scripts to be more consistent and better organized. 
