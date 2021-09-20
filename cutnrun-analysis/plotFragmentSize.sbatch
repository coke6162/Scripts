#!/bin/bash

## Script for running deeptools 
## Date: 30 November 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/scratch/Users/coke6162/yyyymmdd_CUTnRUN_projectName/bwa/filtered \
## outDir=/scratch/Users/coke6162/yyyymmdd_CUTnRUN_projectName/bwa/filtered/fragment_size \
## sbatch --array 0-13 plotFragmentSize.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=4:00:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J plotFragmentSize
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# Define query files
queries=($(ls ${inDir}/*_filtered.sorted.bam | xargs -n 1 basename | sed 's/_filtered.sorted.bam//g' | grep -v "lessThan120" | grep -v "moreThan120"))

# Load modules
module load singularity

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
numThreads=4

# Get fragment size 
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Starting deeptools bamPEFragmentSize...")

singularity exec --bind /scratch ${deeptools} \
bamPEFragmentSize \
--bamfiles ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--table ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_metrics.txt \
--histogram ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_hist.png \
--samplesLabel ${queries[$SLURM_ARRAY_TASK_ID]} \
--plotTitle "Fragment size of PE sorted bam" \
--outRawFragmentLengths ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_rawFragLengths.tab \
--maxFragmentLength 1000 \
--numberOfProcessors ${numThreads} 

echo $(date +"[%b %d %H:%M:%S] Done!")
