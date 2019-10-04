#!/bin/bash

## Script for extracting repeat coordinates
## Date: 1 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/genome-setup
## genomeName=bosTau9 \
## sbatch eachRep_extract.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=0:15:00
#SBATCH --mem=1G

# Job name and output
#SBATCH -J eachRep_extract
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coke6162@colorado.edu

# Define variables
genomeDir=/Users/CL_Shared/db/genomes/${genomeName}

# Create an intermediate file giving only repeat names
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Creating intermediate file containing only repeat names...")

cat \
${genomeDir}/repeats/${genomeName}_allRep_names.txt \
| awk '{print $1}' \
> ${genomeDir}/repeats/${genomeName}_allRep_names.tmp

# Extract repeat coordinates
echo $(date +"[%b %d %H:%M:%S] Extracting repeat coordinates")

while read -r name
do
	grep ${name} ${genomeDir}/repeats/${genomeName}_repeats.bed > ${genomeDir}/repeats/eachRep/tmp/${name}.bed
done < ${genomeDir}/repeats/${genomeName}_allRep_names.tmp

echo $(date +"[%b %d %H:%M:%S] Done")

# Delete temporary file
echo $(date +"[%b %d %H:%M:%S] Removing temporary file")

rm ${genomeDir}/repeats/${genomeName}_allRep_names.tmp

echo $(date +"[%b %d %H:%M:%S] Done")
