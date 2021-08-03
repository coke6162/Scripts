#!/bin/bash

## Script for running MELT Preprocess
## Date: 3 August 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## inDir=/scratch/Users/identikey/yyyymmdd_projectName/bwa/picard/marked \
## outDir=/scratch/Users/identikey/yyyymmdd_projectName/melt \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.fa \
## sbatch --array 0-3 MELT_1_Preprocess.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J Preprocess
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=identikey@colorado.edu

# Set constant variables
numThreads=1

# Define path variables
tmpDir=${inDir}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(ls ${inDir}/*.bam | xargs -n 1 basename))

# Load modules
module load bowtie
module load samtools

# Make temporary directory
pwd; hostname; date

echo "MELT Preprocess version: "$(java -jar /Shares/CL_Shared/programs/MELTv2.2.2/MELT.jar Preprocess --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directories..."

mkdir -p ${tmpDir}

# Run MELT
echo $(date +"[%b %d %H:%M:%S] Preprocessing bams with MELT...")

java \
-Xmx4G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /Shares/CL_Shared/programs/MELTv2.2.2/MELT.jar \
Preprocess \
-bamfile ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-h ${fastaFile}

echo $(date +"[%b %d %H:%M:%S] Done")