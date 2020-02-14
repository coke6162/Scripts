#!/bin/bash

## Script for running MELT Preprocess
## Date: 11 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=20200204_wgs_bosTau9_USMARC_BL3_hybrid \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.main.fa \
## sbatch --array 0-6 MELT_1_Preprocess.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J Preprocess
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=8
memPerThread=4G

# Define path variables
inDir=${parentDir}/${projectName}/gatk/recalibrated
outDir=${parentDir}/${projectName}/melt
tmpDir=${parentDir}/${projectName}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(ls ${inDir}/*_merged_marked_realigned_recalibrated.bam | xargs -n 1 basename))

# Load modules
module load bowtie
module load samtools

# Make temporary directory
pwd; hostname; date

echo "MELT Preprocess version: "$(java -jar /Shares/CL_Shared/programs/MELTv2.1.5/MELT.jar Preprocess --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directories..."

mkdir -p ${tmpDir}
mkdir -p ${outDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Sort with samtools
echo $(date +"[%b %d %H:%M:%S] Sorting recalibrated bams...")

samtools sort \
-@ ${numThreads} \
-T ${outDir}/samtools_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} \
-m ${memPerThread} \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}.sorted.bam

# Index with samtools
echo $(date +"[%b %d %H:%M:%S] Indexing sorted, recalibrated bams...")

samtools index \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}.sorted.bam.bai

# Run MELT
echo $(date +"[%b %d %H:%M:%S] Preprocessing bams with MELT...")

java \
-Xmx32G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /Shares/CL_Shared/programs/MELTv2.1.5/MELT.jar \
Preprocess \
-bamfile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bam}.sorted.bam \
-h ${fastaFile}

echo $(date +"[%b %d %H:%M:%S] Done")