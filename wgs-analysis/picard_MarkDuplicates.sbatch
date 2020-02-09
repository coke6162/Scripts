#!/bin/bash

## Script for running picard MarkDuplicates
## Date: 7 Februrary 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## sbatch --array 0-35 picard_MarkDuplicates.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=72:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J MarkDuplicates
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=16

# Define path variables
inDir=${parentDir}/${projectName}/bwa/picard/merged
outDir=${parentDir}/${projectName}/bwa/picard/marked
tmpDir=${parentDir}/${projectName}/java_tmp

# Define query files
queries=($(ls ${inDir}/*_merged.bam | xargs -n 1 basename))

# Make temporary directory
pwd; hostname; date

echo "Picard MarkDuplicates version: "$(java -jar /opt/picard/2.6.0/picard-2.6.0.jar MarkDuplicates --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}/${SLURM_JOBID}

# Run picard
echo $(date +"[%b %d %H:%M:%S] Marking duplicates with picard...")

java \
-Xmx16G \
-Djava.io.tmpdir=${tmpDir}/${SLURM_JOBID} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/picard/2.6.0/picard-2.6.0.jar \
MarkDuplicates \
INPUT=${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
OUTPUT=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_merged.bam}_merged_marked.bam \
METRICS_FILE=${outDir}/metrics/${queries[$SLURM_ARRAY_TASK_ID]%.bam}_metrics.txt \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=STRICT \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500

echo $(date +"[%b %d %H:%M:%S] Done")