#!/bin/bash

## Script for computing bam alignment count matrix
## Date: 16 March 2019 

## Example usage:
## projectName=project20xx_chipseq \
## bedDir=/Users/coke6162/project20xx_chipseq/macs2/filtered_peaks \
## bed=all_project20xx_chipseq.narrowPeaks \
## outDir=/Users/coke6162/project20xx_chipseq/multicov \
## sbatch computeBamCountMatrix.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bamCountMatrix
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load bedtools

# define key variables
bams=$(cat /Users/coke6162/repos/project_files/"$projectName"/bams.txt)

# run the thing
pwd; hostname; date

echo "Starting bedtools..."
echo $(date +"[%b %d %H:%M:%S] Computing matrix of bam alignment counts...")

# Compute a matrix of bam alignment counts for each region of interest

bedtools multicov -bams ${bams} -bed ${bedDir}/${bed} > ${outDir}/bamCounts_${projectName}_d100.tab

echo $(date +"[%b %d %H:%M:%S] Done")