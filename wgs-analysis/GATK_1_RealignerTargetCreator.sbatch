#!/bin/bash

## Script for running GATK RealignerTargetCreator
## Date: 8 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.main.fa \
## vcfFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/vcf/ARS1.2PlusY_BQSR_v3.main.vcf.gz \
## sbatch GATK_1_RealignerTargetCreator.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH --mem=30GB

# Job name and output
#SBATCH -J RealignerTargetCreator
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=16

# Define path variables
inDir=${parentDir}/${projectName}/bwa/picard/marked
outDir=${parentDir}/${projectName}/gatk/realigned
tmpDir=${parentDir}/${projectName}/java_tmp

# Define variables
bams=$(ls ${inDir}/*_merged_marked.bam | sed 's,^,-I ,g' | tr "\n" " ")

# Make temporary directory
pwd; hostname; date

echo "GATK RealignerTargetCreator version: "$(java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar -T RealignerTargetCreator --version)

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Run GATK
echo $(date +"[%b %d %H:%M:%S] Finding reagions that need to be realigned...")

java \
-Xmx30G \
-Djava.io.tmpdir=${tmpDir}/${SLURM_JOBID} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-nt ${numThreads} \
-R ${fastaFile} \
$(echo ${bams}) \
-rf BadCigar \
-o ${outDir}/target_regions_for_IndelRealigner.intervals \
-known ${vcfFile}

echo $(date +"[%b %d %H:%M:%S] Done")