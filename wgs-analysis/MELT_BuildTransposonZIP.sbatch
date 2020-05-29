#!/bin/bash

## Script for running MELT BuildTransposonZIP
## Date: 29 May 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## zipDir=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/repeats/MELT \
## name=BOVA2 \
## error=10 \
## sbatch MELT_BuildTransposonZIP.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J BuildTransposonZIP
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=1

# Load modules
module load bowtie
module load samtools

# Run MELT
pwd; hostname; date

echo "MELT BuildTransposonZip version: "$(java -jar /Shares/CL_Shared/programs/MELTv2.1.5/MELT.jar BuildTransposonZIP --version)

echo "Preparing MEI transposon reference file..."

java \
-Xmx1G \
-jar /Shares/CL_Shared/programs/MELTv2.1.5/MELT.jar \
BuildTransposonZIP \
${zipDir}/${name}.fa \
${zipDir}/${name}.bed \
${name} \
${error}

mv ${name}_MELT.zip ${zipDir}

echo $(date +"[%b %d %H:%M:%S] Done")