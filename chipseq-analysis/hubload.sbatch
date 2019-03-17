#!/bin/bash

## Script for adding tracks with hubload.py
## Date: 16 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_chipseq/bigwig \
## project=project20xx_chipseq \
## trackDb=/Users/coke6162/hub/mm10/trackDb.txt \
## sbatch hubload.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

# Job name and output
#SBATCH -J hubload
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# commands go here
pwd; hostname; date

echo "Project: "${project}
echo "TrackDb: "${trackDb}

echo $(date +"[%b %d %H:%M:%S] Adding tracks to hub...")

for i in ${inDir}/*.bw 
do
	echo "Processing file: "${i}
	/opt/python/2.7.14/bin/python2.7 /scratch/Shares/chuong/software/scripts/hubload.py --input ${i} --supertrack ${project} --trackDb ${trackDb}
	echo $(date +"[%b %d %H:%M:%S] Done")
done

echo $(date +"[%b %d %H:%M:%S] All done")