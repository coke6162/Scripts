#!/bin/bash

## Script for filtering peaks by excluding blacklisted genome regions 
## Date: 16 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/macs2 \
## outDir=/Users/coke6162/project20xx_chipseq/macs2/filtered_peaks \
## sbatch --array 0-42 filterPeaks.sbatch 

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=1:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J filterPeaks
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.narrowPeak | xargs -n 1 basename))

# Define key variables
blacklist=/Users/coke6162/yoshida2019_atacseq/genrich/mm10_blacklisted_regions_ENCFF547MET.bed

# Other options:
# /Users/CL_Shared/db/genomes/hg38/blacklist/hg38.blacklist.bed

# Load modules
module load bedtools

# Remove blacklisted regions with bedtools intersect
pwd; hostname; date

echo "bedtools version: "$(bedtools --version)
echo "Processing peak file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Starting bedtools intersect..."

bedtools intersect -v -a ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} -b ${blacklist} > ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.filtered

echo $(date +"[%b %d %H:%M:%S] Done")