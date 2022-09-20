#!/bin/bash

## Script for running deeptools 
## Date: 10 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/3_bams/control \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/reports/bam_fragment_size \
## sbatch --array 0-0 plot_fragment_size.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --time=4:00:00
#SBATCH --mem=2GB

# Job name and output
#SBATCH -J plotFragmentSize
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*_filtered.sorted.bam | xargs -n 1 basename | sed 's/_filtered.sorted.bam//g'))

# load modules
module purge
module load python/3.6.3

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img
numThreads=2

# Define temporary directory
tmpDir=${outDir}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directories
pwd; hostname; date

echo "deeptools version: "$(deeptools --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directories..."

mkdir -p ${tmpDir}/deeptools

# Define temporary directory for deeptools
# Note that if you do not set this deeptools will write temporary files to /var and your jobs will likely be killed
TMPDIR="${tmpDir}/deeptools"
export TMPDIR

# Get fragment size 
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Starting deeptools bamPEFragmentSize...")

bamPEFragmentSize \
--bamfiles ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--table ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_metrics.txt \
--histogram ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_hist.png \
--samplesLabel ${queries[$SLURM_ARRAY_TASK_ID]} \
--plotTitle "Fragment size of PE sorted bam" \
--outRawFragmentLengths ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_rawFragLengths.tab \
--maxFragmentLength 1000 \
--numberOfProcessors ${numThreads} 

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing intermediate files...")

rm -rf ${tmpDir}/deeptools
rm -rf ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")