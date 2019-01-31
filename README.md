# scripts
Scripts used for re-analysis of published datasets. 

Over time, (hopefully) this repos will contain generic scripts for the following dataset types:
* ChIP-seq
* RNA-seq
* ATAC-seq
* PRO-seq

These scripts are designed to work on CU Boulder's native computing system, Fiji, which uses Slurm as a resource manager and job scheduler. 

These scripts are intentionally generic - each should contain constant variables that, ideally, should not need to be changed between projects. Project specific path variables, such as input and output directories, will be specified at the command line when submitting the script.