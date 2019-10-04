#!/bin/bash

## Script for uploading bigwigs to UCSC genome browser
## Date: 13 September 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/project20xx_cutnrun/bigwig \
## orderDir=/Users/coke6162/repos/project_files/project20xx_cutnrun \
## projectName=project20xx_cutnrun \
## trackDb=/Users/coke6162/hub/bosTau9/trackDb.txt \
## sbatch hubload.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=4GB

# Job name and output
#SBATCH -J hubload
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Run hubload
pwd; hostname; date

echo $(date +"[%b %d %H:%M:%S] Running hubload...")

for i in ${inDir}/*.bw
do
echo "Processing file: "${i}
/opt/python/2.7.14/bin/python2.7 /Shares/CL_Shared/scripts/ed/hubload.py \
--input ${i} \
--supertrack ${projectName} \
--trackDb ${trackDb}
echo $(date +"[%b %d %H:%M:%S] Done")
done

for i in ${inDir}/*.sorted.bam
do
echo "Processing file: "${i}
/opt/python/2.7.14/bin/python2.7 /Shares/CL_Shared/scripts/ed/hubload.py \
--input ${i} \
--supertrack ${projectName} \
--trackDb ${trackDb}
echo $(date +"[%b %d %H:%M:%S] Done")
done

echo $(date +"[%b %d %H:%M:%S] Done")