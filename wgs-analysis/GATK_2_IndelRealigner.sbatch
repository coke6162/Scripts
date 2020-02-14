#!/bin/bash

## Script for running GATK IndelRealigner
## Date: 8 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.main.fa \
## vcfFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/vcf/ARS1.2PlusY_BQSR_v3.main.vcf.gz \
## sbatch --array 0-6 GATK_2_IndelRealigner.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=48:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J IndelRealigner
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=1

# Define path variables
inDir=${parentDir}/${projectName}/bwa/picard/marked
outDir=${parentDir}/${projectName}/gatk/realigned
tmpDir=${parentDir}/${projectName}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(ls ${inDir}/*_merged_marked.bam | xargs -n 1 basename))

# Make temporary directory
pwd; hostname; date

echo "GATK IndelRealigner version: "$(java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar -T IndelRealigner --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Run GATK
echo $(date +"[%b %d %H:%M:%S] Finding regions that need to be realigned...")

java \
-Xmx16G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ${fastaFile} \
-targetIntervals ${outDir}/target_regions_for_IndelRealigner.intervals \
-I ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged_marked.bam}_merged_marked_realigned.bam \
-known ${vcfFile} \
-rf BadCigar

echo $(date +"[%b %d %H:%M:%S] Done")