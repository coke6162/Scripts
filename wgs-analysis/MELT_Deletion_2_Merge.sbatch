#!/bin/bash

## Script for running MELT Deletion-Merge
## Date: 24 February 2020 

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=20200204_wgs_bosTau9_USMARC_BL3_hybrid \
## zipDir=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/repeats/MELT \
## fastaFile=/Shares/CL_Shared/db/genomes/bosTau9/1000_Bull_Genomes_Project/fa/ARS-UCD1.2_Btau5.0.1Y.main.fa \
## meiList=/Users/coke6162/repos/project_files/20200204_wgs_bosTau9_USMARC_BL3_hybrid/mei_list.txt \
## sbatch --array 0-5 MELT_Deletion_2_Merge.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=14-0:00:00
#SBATCH --mem=30GB

# Job name and output
#SBATCH -J Deletion-Merge
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
inDir=${parentDir}/${projectName}/melt/deletion
outDir=${parentDir}/${projectName}/melt/deletion
tmpDir=${parentDir}/${projectName}/java_tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Define query files
queries=($(cat ${meiList} | xargs -n 1 basename | sed 's/_MELT.zip//g'))

# Make directories
pwd; hostname; date

echo "Making temporary directory..."

mkdir -p ${tmpDir}

# Run MELT
echo $(date +"[%b %d %H:%M:%S] Running MELT Deletion (Step 2: Merge)...")

java \
-Xmx30G \
-Djava.io.tmpdir=${tmpDir} \
-XX:ParallelGCThreads=${numThreads} \
-jar /Shares/CL_Shared/programs/MELTv2.1.5/MELT.jar \
Deletion-Merge \
-mergelist ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}/list.of.outputs.txt \
-bed ${zipDir}/${queries[$SLURM_ARRAY_TASK_ID]}.bed \
-h ${fastaFile} \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done")