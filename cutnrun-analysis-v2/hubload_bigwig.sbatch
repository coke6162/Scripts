#!/bin/bash

## Script for adding bigwig tracks with hubload.py
## Date: 11 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/7_bigwigs \
## trackName=yyyymmdd_projectName \
## trackdb=/Users/coke6162/hub/bosTau9/trackDb.txt \
## sbatch hubload_bigwig.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J hubload_bigwig
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load modules
module load python/2.7.14

# Add bigwigs to hub
pwd; hostname; date

echo "Track name: "${trackName}
echo "TrackDb: "${trackdb}

echo $(date +"[%b %d %H:%M:%S] Adding bigwig tracks to hub...")

for i in ${inDir}/*.bw 
do
	echo "Processing file: "${i}
	hubload.py --input ${i} --supertrack ${trackName} --trackDb ${trackdb} --autoScale off
	echo $(date +"[%b %d %H:%M:%S] Done")
done

echo $(date +"[%b %d %H:%M:%S] Done")