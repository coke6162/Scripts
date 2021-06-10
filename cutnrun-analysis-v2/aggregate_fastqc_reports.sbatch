#!/bin/bash

## Script for running multiqc
## Date: 10 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/reports/raw_fastqc \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/reports/raw_multiqc \
## sbatch aggregate_fastqc_reports.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=0:30:00
#SBATCH --mem=1GB

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
# Note that this assumes you are running on /scratch - if not, change '--bind' to reference your filesystem
pwd; hostname; date

echo "Singularity version: "$(singularity --version)
echo "MultiQC version: "$(singularity exec --bind /scratch ${multiqc} multiqc --version)
echo $(date +"[%b %d %H:%M:%S] Running MultiQC...")

singularity exec --bind /scratch ${multiqc} \
multiqc \
${inDir}/*.zip \
--outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")