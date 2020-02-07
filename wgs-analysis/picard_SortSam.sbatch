#!/bin/bash

## Script for running picard SortSam
## Date: 6 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## projectName=yyyymmdd_wgs_projectName \
## sbatch --array 0-35 picard_SortSam.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J SortSam
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8

# Define path variables
inDir=$HOME/${projectName}/bwa/picard/sorted
outDir=$HOME/${projectName}/bwa/picard/merged
tmpDir=$HOME/${projectName}/java_tmp

# Define query files
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename))

# Make temporary directory
pwd; hostname; date

echo "Picard SortSam version: "$(java -jar /opt/picard/2.6.0/picard-2.6.0.jar SortSam --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Run picard
echo $(date +"[%b %d %H:%M:%S] Sorting bam files with picard...")

java \
-Xmx8G \
-Djava.io.tmpdir=${tmpDir}/${SLURM_JOBID} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/picard/2.6.0/picard-2.6.0.jar \
SortSam \
INPUT=${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
OUTPUT=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=STRICT

echo $(date +"[%b %d %H:%M:%S] Done")