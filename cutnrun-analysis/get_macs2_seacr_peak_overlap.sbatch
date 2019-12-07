#!/bin/bash

## Script for getting the intersection between macs2 and seacr peaks
## Date: December 5 2019

## Example usage:
## cd /Users/coke6162/repos/project_files/20191111_CUTnRUN_MDCK
## macsDir=/Users/coke6162/20191111_CUTnRUN_MDCK/macs2 \
## seacrDir=/Users/coke6162/20191111_CUTnRUN_MDCK/seacr/stringent/slopSummits/b75 \
## outDir=/Users/coke6162/20191111_CUTnRUN_MDCK/macs2/seacrOverlap/macs2_summits_seacr_150bp \
## bash run_get_macs2_seacr_peak_overlap.sh

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=8G

# Job name and output
#SBATCH -J overlap_macs2_seacr
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# load modules
module load bedtools

# run the thing
pwd; hostname; date

echo "Identifying macs2 summits that overlap SEACR maximum signal regions..."
echo "Bedtools version: "$(bedtools --version)

echo "Processing macs2 summits file: "${macsDir}/${macsFile}
echo "Processing seacr maxSignalRegion file: "${seacrDir}/${seacrFile}

bedtools intersect \
-a ${macsDir}/${macsFile} \
-b ${seacrDir}/${seacrFile} \
-u \
> ${outDir}/${outFile}

echo $(date +"[%b %d %H:%M:%S] Done")