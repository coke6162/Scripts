#!/bin/bash

## Script for running bwa index
## Date: 17 September 2019 

## Example usage:
## cd /Users/coke6162/repos/scripts/genome-setup
## inDir=/Shares/CL_Shared/db/genomes/bosTau9/fa \
## outDir=/Shares/CL_Shared/db/genomes/bosTau9/index/bwa \
## scriptDir=/Users/coke6162/repos/scripts/genome-setup \
## prefix=bosTau9.main \
## fastaFile=bosTau9.main.fa \
## sbatch bwa_index.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=4:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bwa
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Load bwa 
module load bwa

# Run bwa
pwd; hostname; date

echo "Preparing bwa indices..."
echo "bwa version: "$(bwa --version)

bwa index \
-p ${prefix} \
${inDir}/${fastaFile} \

# Move output
pwd; hostname; date

echo "Moving output files..."

mv ${scriptDir}/${prefix}* ${outDir}/

echo $(date +"[%b %d %H:%M:%S] Done")