#!/bin/bash

## Script for running macs2
## Date: September 13 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/project20xx_cutnrun/bwa/filtered \
## outDir=/Users/coke6162/project20xx_cutnrun/macs2 \
## genomeSize=2.66e9 \
## sbatch --array=0-3 macs2_PE.sbatch

## genome sizes not hardcoded: bosTau8 (2.66e9), canFam3 (2.33e09)

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

# Job name and output
#SBATCH -J macs2
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Module load macs2
module load python/2.7.14/MACS/2.1.1

# Define query files
queries=($(ls ${inDir}/*_filtered.sorted.bam | xargs -n 1 basename | sed 's/_filtered.sorted.bam//g' | uniq))

# Call peaks with macs2
pwd; hostname; date

echo "Calling peaks with macs2"
echo "macs2 version: "$(macs2 --version)

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

macs2 callpeak \
--treatment ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--name ${queries[$SLURM_ARRAY_TASK_ID]} \
--outdir ${outDir} \
--gsize ${genomeSize} \
--SPMR -B --keep-dup all --format BAMPE --call-summits

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '--treatment <file>' - input file
# '--name <name>' - file name; will be used to name output
# '-B' - whether or not to save extended fragment pileup and local lambda tracks at every bp into a bedGraph file (default = false)
# '--SPMR' - save signal per million reads for fragment pileup profiles; requires -B to be set (default = false)
# '--gsize <int>' - effective genome size; can be 1.0e+9 or 1000000000; can also use shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
# '--format <AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE>' - format of tag file (default = AUTO); use '-f BAMPE' if using paired end data
# '--outdir <outDir>' - output will be written here
# '--keep-dup' - specifies whether macs2 will retain or discard duplicates