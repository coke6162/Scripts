#!/bin/bash

## Script for running featureCounts
## Date: 4 May 2019 

## Example usage:
## cd /Users/coke6162/repos/scripts/rnaseq-analysis
## inDir=/Users/coke6162/project20xx_rnaseq/hisat2 \
## outDir=/Users/coke6162/project20xx_rnaseq/featureCounts \
## projectDir=/Users/coke6162/repos/project_files/project20xx_rnaseq \
## annotationFile=/Users/CL_Shared/db/genomes/bosTau8/annotation/Bos_taurus.UMD3.1.94.gtf \
## sbatch featureCounts.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J featureCounts
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define constant variables
numThreads=8

# Load subread
module load subread

# Run featureCounts, outputting a single table containing count information for all samples
# Note: this does not include any normalization 
pwd; hostname; date

echo "featureCounts version: "$(featureCounts -v)
echo "Starting featureCounts gene-level summarization..."

featureCounts \
-p \
-O \
-T ${numThreads} \
-s 0 \
-t exon \
-g gene_id \
-a ${annotationFile} \
-o ${outDir}/featureCounts_s0.txt \
$(cat ${projectDir}/bam_order.txt)

featureCounts \
-p \
-O \
-T ${numThreads} \
-s 1 \
-t exon \
-g gene_id \
-a ${annotationFile} \
-o ${outDir}/featureCounts_s1.txt \
$(cat ${projectDir}/bam_order.txt)

featureCounts \
-p \
-O \
-T ${numThreads} \
-s 2 \
-t exon \
-g gene_id \
-a ${annotationFile} \
-o ${outDir}/featureCounts_s2.txt \
$(cat ${projectDir}/bam_order.txt)

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-O' - assign reads to all their overlapping meta-features
# '-T <int>' - number of threads; 1 by default
# '-s <int>' - perform strand-specific read counting; options: 0 (unstranded), 1 (stranded), and 2 (reversely stranded); 0 by default
# '-g <option>' - specify attribute type in GTF annotation; 'gene_id' by default
# '-a <file.gtf' - name of annotation file; gtf format by default