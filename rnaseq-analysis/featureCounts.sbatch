#!/bin/bash

## Script for running featureCounts
## Date: 30 March 2019 

## Example usage:
## inDir=/Users/coke6162/project20xx_rnaseq/hisat2 \
## outDir=/Users/coke6162/project20xx_rnaseq/featureCounts \
## annotationFile=/scratch/Shares/chuong/genomes/mm10/gencode.vM18.annotation.gtf
## sbatch featureCounts.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J featureCounts
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%j.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%j.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Load subread
module load subread

# Run featureCounts, outputting a single table containing count information for all samples
# Note: this does not include any normalization 
pwd; hostname; date

echo "featureCounts version: "$(featureCounts -v)
echo "Starting featureCounts gene-level summarization..."

for i in {0..2} # "-s 0" for unstranded, "-s 1" for stranded, "-s 2" for reversely stranded
do
featureCounts -O -T 1 -s ${i} -t exon -g gene_id -a ${annotationFile} -o ${outDir}/featureCounts_s${i}.txt ${inDir}/*.sorted.uniq.bam &
done
wait

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '-O' - assign reads to all their overlapping meta-features
# '-T <int>' - number of threads; 1 by default
# '-s <int>' - perform strand-specific read counting; options: 0 (unstranded), 1 (stranded), and 2 (reversely stranded); 0 by default
# '-g <option>' - specify attribute type in GTF annotation; 'gene_id' by default
# '-a <file.gtf' - name of annotation file; gtf format by default