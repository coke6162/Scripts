#!/bin/bash

## Script for running MACS2
## Date: 10 July 2019 

## Example usage:
## cd /Users/coke6162/repos/project_files/project20xx_chipseq
## inDir=/Users/coke6162/project20xx_chipseq/bwa \
## bash runMacs2.sh

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J macs2_callpeak
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load MACS2
module load python/2.7.14/MACS/2.1.1

# Run macs2
pwd; hostname; date

echo "MACS2 version: "$(macs2 --version)
echo "Processing treatment file: "${treatFile}
echo "Control file: "${controlFile}
echo "Genome size: "${genomeSize}
echo $(date +"[%b %d %H:%M:%S] Starting MACS2 peak calling...")

macs2 callpeak \
-t ${treatFile} \
-c ${controlFile} \
-n ${name} \
-B --SPMR \
-g ${genomeSize} \
-f BAM \
--outdir ${outDir}

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-t <file>' - ChIP-seq treatment file
# '-c <file>' - ChIP-seq input file
# '-n <name>' - file name; will be used to name output
# '-B' - whether or not to save extended fragment pileup and local lambda tracks at every bp into a bedGraph file (default = false)
# '--SPMR' - save signal per million reads for fragment pileup profiles; requires -B to be set (default = false)
# '-g <int>' - effective genome size; can be 1.0e+9 or 1000000000; can also use shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
# '-f <AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE>' - format of tag file (default = AUTO); use '-f BAMPE' if using paired end data
# '--outdir <outDir>' - output will be written here