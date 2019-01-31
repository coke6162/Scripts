#!/bin/bash
#SBATCH -p short                   # Partition or queue
#SBATCH --job-name=project20xx_rnaseq_hisat2         # Job name
#SBATCH --mail-type=END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=coke6162@colorado.edu     # Where to send mail
#SBATCH --ntasks=1                 # Run on a single CPU
#SBATCH --cpus-per-task=8         # cpus
#SBATCH --mem=16gb                # Memory limit
#SBATCH --time=3:00:00            # Time limit hrs:min:sec
#SBATCH --output=/Users/coke6162/sbatch_script/output/project20xx_rnaseq_hisat2_%A_%a.out      # Standard output log
#SBATCH --error=/Users/coke6162/sbatch_script/error/project20xx_rnaseq_hisat2_%A_%a.err       # Standard error log

# Set constant variables
refGenome=mm10 # mm10, hg38, canFam3, bosTau8
indexFileName=mm10
rnaStrandness=FR # specify whether the library is forward (FR) or reverse stranded (RF)
numThreads=8

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls ${inDir}/*_paired.fastq.gz | xargs -n 1 basename | sed --expression='s/_R1_paired.fastq.gz//g' \
	| sed --expression='s/_R2_paired.fastq.gz//g' | sed --expression='s/_R1_unpaired.fastq.gz//g' \
	| sed --expression='s/_R2_unpaired.fastq.gz//g' | uniq))

# Load HISAT2 and samtools
module load hisat2 samtools

# Run HISAT2
pwd; hostname; date

echo "Starting HISAT2 alignment..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

hisat2 -p ${numThreads} --rna-strandness ${rnaStrandness} --no-softclip -x /scratch/Shares/chuong/genomes/${refGenome}/${indexFileName} \
-1 ${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz -2 ${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_paired.fastq.gz}.sorted.uniq.bam
wait

# Create .sorted.uniq.bai index files for upload to UCSC browser
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_paired.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

# Command line entry
### inDir=/Users/coke6162/project20xx_rnaseq/trimmomatic \
### outDir=/Users/coke6162/project20xx_rnaseq/hisat2 \
### sbatch --array 0-12 /Users/coke6162/project20xx_rnaseq/scripts/hisat2.sbatch
