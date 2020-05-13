#!/bin/bash

## Script for running picard MergeSamFiles
## Date: 7 April 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=20200204_wgs_bosTau9_USMARC_BL3_hybrid \
## sbatch --array 0-79 picard_MergeSamFiles.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J MergeSamFiles
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Set constant variables
numThreads=16

# Define path variables
inDir=${parentDir}/${projectName}/bwa/picard/sorted
outDir=${parentDir}/${projectName}/bwa/picard/merged
tmpDir=${parentDir}/${projectName}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query samples
# This assumes your samples are uniquely identified at the start, followed by an underscore
# i.e. "Limousin1_run1_lane1.sorted.bam"
queries=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename | cut -d_ -f1 | sort | uniq))

# Make temporary directory
pwd; hostname; date

echo "Picard MergeSamFiles version: "$(java -jar /opt/picard/2.6.0/picard-2.6.0.jar MergeSamFiles --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directory..."

mkdir -p ${tmpDir}

# Create a text file specifying all bam files to merge
# This will put "INPUT=" in front of each file (needed to run Picard MergeSamFiles)
# This works if all your samples share a common suffix & end in ".bam"
echo $(date +"[%b %d %H:%M:%S] Creating a txt file giving names of all files to merge...")

ls ${inDir} | grep "^${queries[$SLURM_ARRAY_TASK_ID]}.*bam$" - > ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt.tmp
sed 's,^,INPUT='"${inDir}"'\/,g' ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt.tmp > ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt
rm ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt.tmp

# Run Picard MergeSamFiles
echo $(date +"[%b %d %H:%M:%S] Merging bam files with picard...")

java \
-Xmx8G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /opt/picard/2.6.0/picard-2.6.0.jar \
MergeSamFiles \
$(cat ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merge_list.txt) \
OUTPUT=${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_merged.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
USE_THREADING=true \
MERGE_SEQUENCE_DICTIONARIES=true \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=STRICT

echo $(date +"[%b %d %H:%M:%S] Done")