#!/bin/bash

## Script for running picard MergeSamFiles
## Date: 6 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## sample=Limousin2 \
## sbatch picard_MergeSamFiles.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J MergeSamFiles
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=16

# Define path variables
inDir=${parentDir}/${projectName}/bwa/picard/sorted
outDir=${parentDir}/${projectName}/bwa/picard/merged
tmpDir=${parentDir}/${projectName}/java_tmp

# Make temporary directory
pwd; hostname; date

echo "Picard MergeSamFiles version: "$(java -jar /opt/picard/2.6.0/picard-2.6.0.jar MergeSamFiles --version)

echo "Processing sample: "${inDir}/${sample}

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Create a text file specifying all bam files to merge
# This will put "INPUT=" in front of each file (needed to run Picard MergeSamFiles)
# This works if all your samples share a common suffix & end in ".bam"
echo $(date +"[%b %d %H:%M:%S] Creating a txt file giving names of all files to merge...")

ls ${inDir} | grep "^${sample}.*bam$" - > ${inDir}/${sample}_merge_list.txt.tmp
sed 's,^,INPUT='"${inDir}"'\/,g' ${inDir}/${sample}_merge_list.txt.tmp > ${inDir}/${sample}_merge_list.txt
rm ${inDir}/${sample}_merge_list.txt.tmp

# Run Picard MergeSamFiles
echo $(date +"[%b %d %H:%M:%S] Merging bam files with picard...")

java \
-Xmx8G \
-Djava.io.tmpdir=${tmpDir}/${SLURM_JOBID} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/picard/2.6.0/picard-2.6.0.jar \
MergeSamFiles \
$(cat ${inDir}/${sample}_merge_list.txt) \
OUTPUT=${outDir}/${sample}_merged.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
USE_THREADING=true \
MERGE_SEQUENCE_DICTIONARIES=true \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=STRICT

echo $(date +"[%b %d %H:%M:%S] Done")