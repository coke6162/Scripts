#!/bin/bash

## Script for subsetting bam files for fragments <= 150bp
## Date: 10 June 2021

## Example usage:
## cd /Users/coke6162/repos/scripts/cutnrun-analysis-v2
## inDir=/scratch/Users/coke6162/yyyymmdd_projectName/3_bams/tf \
## outDir=/scratch/Users/coke6162/yyyymmdd_projectName/3_bams/tf \
## reportsDir=/scratch/Users/coke6162/yyyymmdd_projectName/reports/alignment_stats \
## sbatch --array 0-0 subset_bam_by_fragment_length_150bp.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=12:00:00
#SBATCH --mem=6GB

# Job name and output
#SBATCH -J subset_bam
#SBATCH -o /Users/%u/sbatch_script/output/slurm-%A_%a.out
#SBATCH -e /Users/%u/sbatch_script/error/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=coke6162@colorado.edu

# Define query files
queries=($(ls ${inDir}/*_filtered.sorted.bam | xargs -n 1 basename | grep -v "lessThan150" | grep -v "greaterThan151" | sed 's/_filtered.sorted.bam//g' | uniq))

# Load modules
module purge
module load samtools python/3.6.3

# Define key variables
numThreads=2

# Define temporary directory
tmpDir=${outDir}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Make temporary directories
pwd; hostname; date

echo "deeptools version: "$(deeptools --version)
echo "samtools version: "$(samtools --version)

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo "Making temporary directories..."

mkdir -p ${tmpDir}/samtools
mkdir -p ${tmpDir}/deeptools

# Define temporary directory for deeptools
# Note that if you do not set this deeptools will write temporary files to /var and your jobs will likely be killed
TMPDIR="${tmpDir}/deeptools"
export TMPDIR

## Use "alignmentSieve" to filter reads by size (<= 150bp)
echo $(date +"[%b %d %H:%M:%S] Starting deeptools alignmentSieve...")

echo "Processing sample: "${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam

alignmentSieve \
--bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--maxFragmentLength 150 \
--outFile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.bam.tmp \
--numberOfProcessors ${numThreads}

# Use "alignmentSieve" to filter reads by size (>150bp)
echo $(date +"[%b %d %H:%M:%S] Starting deeptools alignmentSieve...")

alignmentSieve \
--bam ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_filtered.sorted.bam \
--minFragmentLength 151 \
--outFile ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.bam.tmp \
--numberOfProcessors ${numThreads}

# Sort subsetted bams
echo $(date +"[%b %d %H:%M:%S] Sorting subsetted bams...")

samtools sort -@ ${numThreads} -T ${tmpDir} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.bam.tmp \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.sorted.bam

samtools sort -@ ${numThreads} -T ${tmpDir} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.bam.tmp \
-o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.sorted.bam

# Index subsetted bams
echo $(date +"[%b %d %H:%M:%S] Indexing subsetted bams...")

samtools index -@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.sorted.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.sorted.bam.bai

samtools index -@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.sorted.bam \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.sorted.bam.bai

# Generate flagstat files from subsetted bam files
echo $(date +"[%b %d %H:%M:%S] Making flagstat file from subsetted bams...")

samtools flagstat \
-@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.sorted.bam \
> ${reportsDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered_flagstat.txt

samtools flagstat -@ ${numThreads} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.sorted.bam \
> ${reportsDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered_flagstat.txt

# Remove intermediate files
echo $(date +"[%b %d %H:%M:%S] Removing intermediate files...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_lessThan150_filtered.bam.tmp
rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_greaterThan151_filtered.bam.tmp
rm -rf ${tmpDir}/samtools
rm -rf ${tmpDir}/deeptools
rm -r ${tmpDir}

echo $(date +"[%b %d %H:%M:%S] Done")