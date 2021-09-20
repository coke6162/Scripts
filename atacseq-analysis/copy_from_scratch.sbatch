#!/bin/bash

## Script to copy data from /scratch
## Date: 4 May 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## scratch_project=/scratch/Users/coke6162/20210126_ATAC-seq_Johnston2021 \
## sbatch copy_from_scratch.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=2:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J cp
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# run the thing
pwd; hostname; date

cp -r ${scratch_project} /Users/coke6162

echo $(date +"[%b %d %H:%M:%S] Done")