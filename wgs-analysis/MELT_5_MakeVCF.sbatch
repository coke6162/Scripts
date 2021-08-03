#!/bin/bash

## Script for running MELT MakeVCF (Step 5)
## Date: 3 August 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## meltDir=/scratch/Users/coke6162/yyyymmdd_projectName/melt \
## zipDir=/scratch/Users/coke6162/yyyymmdd_projectName/melt/zip \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.fa \
## meiList=/scratch/Users/coke6162/yyyymmdd_projectName/melt/mei_list.txt \
## sbatch --array 0-5 MELT_4_MakeVCF.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J MakeVCF
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=1

# Load modules
module load bowtie
module load samtools

# Define path variables
tmpDir=${meltDir}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(cat ${meiList} | xargs -n 1 basename | sed 's/_MELT.zip//g'))

# Make directories
pwd; hostname; date

echo "Making temporary directory..."

mkdir -p ${tmpDir}

# Run MELT
echo $(date +"[%b %d %H:%M:%S] Running MELT MakeVCF...")

java \
-Xmx32G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /Shares/CL_Shared/programs/MELTv2.2.2/MELT.jar \
MakeVCF \
-genotypingdir ${meltDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-p ${meltDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-o ${meltDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-w ${meltDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-t ${zipDir}/${queries[$SLURM_ARRAY_TASK_ID]}_MELT.zip  \
-h ${fastaFile}

echo $(date +"[%b %d %H:%M:%S] Done")