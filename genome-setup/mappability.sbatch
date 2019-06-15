#!/bin/bash

## Script for making mappability tracks with gem
## Date: 22 May 2019

## Example usage:
## cd /Users/coke6162/scripts/genome-setup
## inputFasta=/Users/CL_Shared/db/genomes/bosTau9/fa/bosTau9.fa \
## chromSizes=/Users/CL_Shared/db/genomes/bosTau9/fa/bosTau9.chrom.sizes \
## name=bosTau8 \
## readLength=100
## outDir=/Users/coke6162/genomes/bosTau9/mappability/100mer \
## sbatch mappability.sbatch

# General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=20-0:00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J gem_mappability
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# define path to gem executable
export PATH="$PATH:/scratch/Shares/chuong/software/gem"

# Define key variables
numCPUs=16

# Create a gem index
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Starting gem indexer...")

gem-indexer \
-i ${inputFasta} \
-o ${outDir}/${name} \
-T ${numCPUs}

# Generate mappability files for different read lengths
echo $(date +"[%b %d %H:%M:%S] Generate mappability files for a read length of ")${readLength}

gem-mappability \
-I ${outDir}/${name}.gem \
-l ${readLength} \
-o ${outDir}/${name}.${readLength}mer \
-T ${numCPUs}

# convert GEM format to wiggle format
echo $(date +"[%b %d %H:%M:%S] convert GEM format to wiggle format...")

gem-2-wig \
-I ${outDir}/${name}.gem \
-i ${outDir}/${name}.${readLength}mer.mappability \
-o ${outDir}/${name}.${readLength}mer

# convert wiggle to bigwig
echo $(date +"[%b %d %H:%M:%S] convert wiggle to bigwig...")

wigToBigWig \
${outDir}/${name}.${readLength}mer.wig \
${chromSizes} \
${outDir}/${name}.${readLength}mer.bw 

echo $(date +"[%b %d %H:%M:%S] Done")