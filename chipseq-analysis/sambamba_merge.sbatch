#!/bin/bash

## Script for running sambamba merge
## Date: 16 March 2019 

## Example usage:
## bash sambamba_merge.sh

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J sambamba_markdup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load samtools
module load samtools

# Define query files
queries=($(ls ${inDir}/*.sorted.uniq.bam | xargs -n 1 basename))

# Merge technical replicates with sambamba merge
# Might need to modify to suit your needs (i.e. two inFiles vs four)
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Running sambamba merge...")

sambamba merge --nthreads=8 ${outDir}/${outFile} ${inDir}/${inFile1} ${inDir}/${inFile2} ${inDir}/${inFile3} ${inDir}/${inFile4}

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-nthreads=<int>' - number of threads to use