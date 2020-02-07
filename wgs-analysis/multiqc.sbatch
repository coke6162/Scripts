#!/bin/bash

## Script for running multiqc
## Date: 6 February 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/wgs-analysis
## parentDir=/scratch/Users/coke6162 \
## projectName=yyyymmdd_wgs_projectName \
## sbatch multiqc.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -J multiqc
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load singularity

# define key variables
multiqc=/scratch/Shares/public/singularity/multiqc-1.7.img

# Define path variables
inDir=${parentDir}/${projectName}/fastq/fastqc
outDir=${parentDir}/${projectName}/fastq/fastqc/multiqc

# Run MultiQC
pwd; hostname; date

echo "Singularity version: "$(singularity --version)
echo "MultiQC version: "$(singularity exec --bind /Users/CL_Shared ${multiqc} multiqc --version)
echo $(date +"[%b %d %H:%M:%S] Running MultiQC...")

singularity exec --bind /Shares/CL_Shared ${multiqc} \
multiqc \
${inDir}/*.zip \
--outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")