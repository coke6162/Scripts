#!/bin/bash

## Script for running MACS2
## Date: 14 March 2019 

## Example usage:
## projectFile=/Users/coke6162/repos/project_files/project20xx_chipseq_macs2.sh
## bash ${projectFile}

## example line from ${projectFile}:
## inDir=/Users/coke6162/project20xx_chipseq/bwa outDir=/Users/coke6162/project20xx_chipseq/macs2 treatFile=stat1_ifng_4h.sorted.uniq.bam controlFile=stat1_input.sorted.uniq.bam -n stat1_ifng_4h -g mm sbatch /Users/coke6162/repos/scripts/chipseq-analysis/macs2_PE.sbatch

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

macs2 callpeak -t ${inDir}/${treatFile} -c ${inDir}/${controlFile} -n ${name} -B --SPMR -g ${genomeSize} -f BAMPE --outdir ${outDir}

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