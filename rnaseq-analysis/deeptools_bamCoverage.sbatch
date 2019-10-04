#!/bin/bash

## Script for running deeptools bamCoverage
## Date: 4 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## inDir=/Users/coke6162/project20xx_rnaseq/hisat2 \
## refDir=/Users/coke6162/repos/project_files/project20xx_rnaseq \
## outDir=/Users/coke6162/project20xx_rnaseq/bigwigs \
## sbatch --array 0-5 deeptools_bamCoverage.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J deeptools_bamCov
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# load modules
module load singularity

# Define key variables
deeptools=/scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img

# Define query files
queries=($(cat ${refDir}/bam_order.txt | xargs -n 1 basename | sed --expression='s/.sorted.uniq.bam//g'))

# Define deeptools parameters
binSize=1
numCPU=16

# Use "bamCoverage" to generate continuous coverage signal from bam files
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo "Starting deeptools..."
echo $(date +"[%b %d %H:%M:%S] Convert bam into bigwig with bamCoverage...")

# Start by generating bigwig from forward strand data
echo $(date +"[%b %d %H:%M:%S] Generate forward strand bigwigs...")

singularity exec --bind /scratch ${deeptools} \
bamCoverage \
-b ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_fwd.bw \
--binSize=${binSize} \
-p ${numCPU} \
--normalizeUsing CPM \
--filterRNAstrand forward \
--ignoreForNormalization chrX chrM \

# Start by generating bigwig from reverse strand data
echo $(date +"[%b %d %H:%M:%S] Generate reverse strand bigwigs...")

singularity exec --bind /scratch ${deeptools} \
bamCoverage \
-b ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.sorted.uniq.bam \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_rev.bw \
--binSize=${binSize} \
-p ${numCPU} \
--normalizeUsing CPM \
--filterRNAstrand reverse \
--ignoreForNormalization chrX chrM

echo $(date +"[%b %d %H:%M:%S] Done!")