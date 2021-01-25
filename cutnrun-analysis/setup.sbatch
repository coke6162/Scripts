#!/bin/bash

## Script to set up directory layout on /scratch
## Date: 25 November 2020

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/scratch/Users/coke6162 \
## projectName=project20xx_CUTnRUN \
## sbatch setup.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=0:10:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J setup
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# run the thing
pwd; hostname; date

echo "Working dir: "${inDir}
echo "Project: "${projectName}
echo $(date +"[%b %d %H:%M:%S] Making project directories...")

# go to the designated work dir
cd ${inDir}

# make a new project dir
mkdir -p ${projectName}

# make necessary subdirs
mkdir -p ${projectName}/fastq
mkdir -p ${projectName}/fastq/fastqc
mkdir -p ${projectName}/fastq/fastqc/multiqc
mkdir -p ${projectName}/bbduk
mkdir -p ${projectName}/bbduk/fastqc
mkdir -p ${projectName}/bbduk/fastqc/multiqc
mkdir -p ${projectName}/bwa
mkdir -p ${projectName}/bwa/flagstat
mkdir -p ${projectName}/bwa/filtered
mkdir -p ${projectName}/bwa/filtered/flagstat
mkdir -p ${projectName}/bwa/filtered/seacrInput
mkdir -p ${projectName}/bwa/filtered/unNormBdg
mkdir -p ${projectName}/macs2
mkdir -p ${projectName}/macs2/FRIP
mkdir -p ${projectName}/seacr
mkdir -p ${projectName}/seacr/relaxed
mkdir -p ${projectName}/seacr/stringent
mkdir -p ${projectName}/seacr/stringent/maximumSignalRegion
mkdir -p ${projectName}/seacr/stringent/slopSummits
mkdir -p ${projectName}/seacr/stringent/slopSummits/b250
mkdir -p ${projectName}/seacr/stringent/slopSummits/b250/fa
mkdir -p ${projectName}/seacr/stringent/slopSummits/b75
mkdir -p ${projectName}/seacr/stringent/slopSummits/b75/fa
mkdir -p ${projectName}/bigwig
mkdir -p ${projectName}/bigwig/macs2
mkdir -p ${projectName}/bigwig/unNormFrag
mkdir -p ${projectName}/repeatAnalysis
mkdir -p ${projectName}/repeatAnalysis/macs2
mkdir -p ${projectName}/repeatAnalysis/seacr
mkdir -p ${projectName}/deeptools
mkdir -p ${projectName}/deeptools/macs2
mkdir -p ${projectName}/deeptools/seacr
mkdir -p ${projectName}/deeptools/seacr/relaxed
mkdir -p ${projectName}/deeptools/seacr/stringent
mkdir -p ${projectName}/giggle
mkdir -p ${projectName}/giggle/macs2
mkdir -p ${projectName}/giggle/macs2/nonDE
mkdir -p ${projectName}/giggle/macs2/nonDE/allRep
mkdir -p ${projectName}/giggle/macs2/nonDE/LTR
mkdir -p ${projectName}/giggle/macs2/nonDE/SINE
mkdir -p ${projectName}/giggle/macs2/nonDE/LINE
mkdir -p ${projectName}/giggle/macs2/nonDE/DNA
mkdir -p ${projectName}/giggle/seacr
mkdir -p ${projectName}/giggle/seacr/stringent
mkdir -p ${projectName}/giggle/seacr/stringent/nonDE
mkdir -p ${projectName}/giggle/seacr/stringent/nonDE/allRep
mkdir -p ${projectName}/giggle/seacr/stringent/nonDE/LTR
mkdir -p ${projectName}/giggle/seacr/stringent/nonDE/SINE
mkdir -p ${projectName}/giggle/seacr/stringent/nonDE/LINE
mkdir -p ${projectName}/giggle/seacr/stringent/nonDE/DNA
mkdir -p ${projectName}/motif_analysis
mkdir -p ${projectName}/motif_analysis/meme-chip
mkdir -p ${projectName}/motif_analysis/meme-chip/seacr
mkdir -p ${projectName}/motif_analysis/meme-chip/seacr/stringent
mkdir -p ${projectName}/motif_analysis/meme-chip/seacr/stringent/slopSummits_b250

echo $(date +"[%b %d %H:%M:%S] Done")