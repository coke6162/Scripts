#!/bin/bash

## Script to copy data from /scratch
## Date: 19 December 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/scratch/Users/coke6162 \
## projectName=project20xx_CUTnRUN \
## destinationDir=/Users/coke6162 \
## sbatch cp.sbatch

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

echo "Working dir: "${inDir}/${projectName}
echo "Copying to: "${destinationDir}

# make necessary directories
echo $(date +"[%b %d %H:%M:%S] Making destination directories...")

mkdir -p ${destinationDir}/${projectName}/fastq/fastqc/multiqc
mkdir -p ${destinationDir}/${projectName}/bwa/filtered/fragment_size
mkdir -p ${destinationDir}/${projectName}/macs2/FRIP
mkdir -p ${destinationDir}/${projectName}/macs2/slopSummits
mkdir -p ${destinationDir}/${projectName}/seacr
mkdir -p ${destinationDir}/${projectName}/bigwig
mkdir -p ${destinationDir}/${projectName}/deeptools
mkdir -p ${destinationDir}/${projectName}/giggle
mkdir -p ${destinationDir}/${projectName}/meme
mkdir -p ${destinationDir}/${projectName}/repeatAnalysis

# go to the designated work dir
echo $(date +"[%b %d %H:%M:%S] Moving to working directory...")

cd ${inDir}/${projectName}

# make necessary subdirs
cp -n -r fastq/ ${destinationDir}/${projectName}/
cp -n bwa/*_chrMreadsFraction.txt ${destinationDir}/${projectName}/bwa
cp -n -r bwa/filtered/fragment_size ${destinationDir}/${projectName}/filtered/fragment_size
cp -n macs2/*_peaks.narrowPeak ${destinationDir}/${projectName}/macs2
cp -n macs2/*_summits.bed ${destinationDir}/${projectName}/macs2
cp -n -r macs2/FRIP ${destinationDir}/${projectName}/macs2
cp -n -r macs2/slopSummits ${destinationDir}/${projectName}/macs2
cp -n -r seacr/ ${destinationDir}/${projectName}/
cp -n -r bigwig/ ${destinationDir}/${projectName}/
cp -n -r deeptools/ ${destinationDir}/${projectName}/
cp -n -r giggle/ ${destinationDir}/${projectName}/
cp -n -r meme/ ${destinationDir}/${projectName}/
cp -n -r repeatAnalysis/ ${destinationDir}/${projectName}/

echo $(date +"[%b %d %H:%M:%S] Done")