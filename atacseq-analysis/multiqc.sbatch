#!/bin/bash

## Script for running multiqc
## Date: 30 April 2019 

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## inDir=/Users/coke6162/project20xx_atacseq/fastqc \
## outDir=/Users/coke6162/project20xx_atacseq/fastqc/multiqc \
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

# Run MultiQC
pwd; hostname; date

echo "MultiQC version: "$(multiqc --version)
echo $(date +"[%b %d %H:%M:%S] Running multiqc...")

multiqc ${inDir}/*.zip --outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")