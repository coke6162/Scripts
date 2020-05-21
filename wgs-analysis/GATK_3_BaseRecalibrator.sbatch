#!/bin/bash

## Script for running GATK BaseRecalibrator
## Date: 8 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.main.fa \
## vcfFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/vcf/ARS1.2PlusY_BQSR_v3.main.vcf.gz \
## sbatch --array 0-6 GATK_3_BaseRecalibrator.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=3-0:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J BaseRecalibrator
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=4

# Define path variables
inDir=${parentDir}/${projectName}/gatk/realigned
outDir=${parentDir}/${projectName}/gatk/recalibrated
tmpDir=${parentDir}/${projectName}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(ls ${inDir}/*_merged_marked_realigned.bam | xargs -n 1 basename))

# Make temporary directory
pwd; hostname; date

echo "GATK BaseRecalibrator version: "$(java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar -T BaseRecalibrator --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Run GATK
echo $(date +"[%b %d %H:%M:%S] Finding reagions that need to be realigned...")

java \
-Xmx16G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-nct ${numThreads} \
-R ${fastaFile} \
-I ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged_marked_realigned.bam}.recaltable \
-knownSites ${vcfFile} \
-bqsrBAQGOP 45 \
-rf BadCigar

echo $(date +"[%b %d %H:%M:%S] Done")