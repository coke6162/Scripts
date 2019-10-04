#!/bin/bash

## Script for adding tracks with hubload.py
## Date: 15 July 2019 

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## inDir=/Users/coke6162/project20xx_atacseq/bigwigs \
## projectName=project20xx_atacseq \
## trackDb=/Users/coke6162/hub/bosTau8/trackDb.txt \
## sbatch hubload.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J hubload
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# commands go here
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Adding tracks to hub...")

echo "Project: "${project}
echo "TrackDb: "${trackDb}

for i in ${inDir}/*.bw 
do
	echo "Processing file: "${i}
	/opt/python/2.7.14/bin/python2.7 /Users/CL_Shared/scripts/ed/hubload.py \
	--input ${i} \
	--supertrack ${projectName} \
	--trackDb ${trackDb}
	echo $(date +"[%b %d %H:%M:%S] Done")
done

echo $(date +"[%b %d %H:%M:%S] Done")