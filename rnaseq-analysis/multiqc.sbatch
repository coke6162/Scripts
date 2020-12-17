#!/bin/bash

## Script for running multiqc
## Date: 28 June 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_rnaseq/fastq/fastqc \
## outDir=/Users/coke6162/project20xx_rnaseq/fastq/fastqc/multiqc \
## sbatch multiqc.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00
#SBATCH --mem=4GB

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

# Run MultiQC
pwd; hostname; date

echo "Singularity version: "$(singularity --version)
echo "MultiQC version: "$(singularity exec --bind /scratch/Users ${multiqc} multiqc --version)
echo $(date +"[%b %d %H:%M:%S] Running MultiQC...")

singularity exec --bind /scratch/Users ${multiqc} \
multiqc \
${inDir}/*.zip \
--outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")
