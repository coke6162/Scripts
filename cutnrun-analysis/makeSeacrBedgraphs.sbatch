#!/bin/bash

## Script for running deeptools 
## Date: 18 November 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis
## inDir=/Users/coke6162/20191111_CUTnRUN_MDCK/bwa/filtered \
## outDir=/Users/coke6162/20191111_CUTnRUN_MDCK/bwa/filtered/seacrInput \
## genomeChrFile=/Shares/CL_Shared/db/genomes/canFam3/fa/canFam3.main.chrom.sizes \
## sbatch --array 0-44 makeSeacrBedgraphs.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J subsetFragSize
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=conor.j.kelly@colorado.edu

# Define query files
queries=($(ls ${inDir}/*.bam | xargs -n 1 basename | sed --expression='s/_filtered.sorted.bam//g' | uniq))

# load modules
module load bedtools samtools

# define key variables
numThreads=8

# Convert to BEDPE & filter out fragments greater than 1kb
pwd; hostname; date

echo "Processing file..."${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam

echo $(date +"[%b %d %H:%M:%S] Converting filtered bam to BEDPE...")

samtools sort \
-@ ${numThreads} \
-n ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted.bam.tmp

samtools fixmate \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted.bam.tmp \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed.bam.tmp

samtools view \
-bf 0x2 ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed.bam.tmp \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed_onlyPaired.bam.tmp

bedtools bamtobed \
-bedpe \
-i ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed_onlyPaired.bam.tmp \
| awk '$1==$4 && $6-$2 < 1000 {print $0}' - \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_d1000.sorted.bed

echo "Sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_stats.txt
echo "Number of reads after sorting by name: "$(samtools view ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted.bam.tmp | wc -l) > ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_stats.txt
echo "Number of reads after mate fixing: "$(samtools view ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed.bam.tmp | wc -l) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_stats.txt
echo "Number of reads after removing improperly paired reads: "$(samtools view ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed_onlyPaired.bam.tmp | wc -l) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_stats.txt
echo "Number of fragments in bed file: "$(cat ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_d1000.sorted.bed | wc -l) >> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_stats.txt

# Select the 5' and 3' coordinates of the read pair to generate a new BED3 file
echo $(date +"[%b %d %H:%M:%S] Extract read pair coordinates to generate BED3...")

cut -f 1,2,6 \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_d1000.sorted.bed \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_d1000_fragments.sorted.bed

# Convert the BED3 files into bedgraphs
echo $(date +"[%b %d %H:%M:%S] Generating bedgraph files...")

sort -k1,1 ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_d1000_fragments.sorted.bed \
| bedtools genomecov \
-bg \
-i - \
-g ${genomeChrFile} \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_d1000_fragments.bdg

# Remove intermediate files
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted.bam.tmp
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed.bam.tmp
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered_resorted_fixed_onlyPaired.bam.tmp

echo $(date +"[%b %d %H:%M:%S] Done")