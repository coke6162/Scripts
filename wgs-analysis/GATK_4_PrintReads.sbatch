#!/bin/bash

## Script for running GATK PrintReads
## Date: 11 Februrary 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.main.fa \
## sbatch --array 0-6 GATK_4_PrintReads.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=48:00:00
#SBATCH --mem=30GB

# Job name and output
#SBATCH -J PrintReads
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=4

# Load java through module
# The default java on Fiji is OpenJDK 11.0.1
# To get GATK/3.7 to work, we need Sun/Oracle java 1.8
module load java

# Define path variables
inDir=${parentDir}/${projectName}/gatk/realigned
outDir=${parentDir}/${projectName}/gatk/recalibrated
tmpDir=${parentDir}/${projectName}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(ls ${inDir}/*_merged_marked_realigned.bam | xargs -n 1 basename | sed 's/_merged_marked_realigned.bam//g'))

# Make temporary directory
pwd; hostname; date

echo "GATK PrintReads version: "$(java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar -T PrintReads --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Run GATK
echo $(date +"[%b %d %H:%M:%S] Printing reads...")

java \
-Xmx30G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
-T PrintReads \
-nct ${numThreads} \
-R ${fastaFile} \
--BQSR ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}.recaltable \
-I ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged_marked_realigned.bam \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged_marked_realigned_recalibrated.bam \
-rf BadCigar

echo $(date +"[%b %d %H:%M:%S] Done")