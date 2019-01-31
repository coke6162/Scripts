#!/bin/bash
#SBATCH -p short                   # Partition or queue
#SBATCH --job-name=project20xx_chipseq_bwa        # Job name
#SBATCH --mail-type=END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=coke6162@colorado.edu     # Where to send mail
#SBATCH --ntasks=1                 # Run on a single CPU
#SBATCH --cpus-per-task=8         # cpus
#SBATCH --mem=64gb                # Memory limit
#SBATCH --time=4:00:00            # Time limit hrs:min:sec
#SBATCH --output=/Users/coke6162/sbatch_script/output/project20xx_chipseq_bwa_%A_%a.out      # Standard output log
#SBATCH --error=/Users/coke6162/sbatch_script/error/project20xx_chipseq_bwa_%A_%a.err       # Standard error log

# Set constant variables
bwaIndexDir=/scratch/Shares/chuong/genomes/hg38 # options: hg38, mm10, bosTau8, canFam3
bwaIndex=hg38.main.fa # options: hg38.main.fa, mm10.fa, bosTau8.fa, canFam3.fa

# Set query files
queries=($(ls ${inDir}/*_trimmed.fastq.gz | xargs -n 1 basename))

# Load bwa and samtools
module load bwa/0.7.15 samtools/1.8

# Run bwa
pwd; hostname; date

echo "Starting bwa alignment..."
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

bwa mem -t 8 ${bwaIndexDir}/${bwaIndex} ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]} \
| samtools view -Sb -q 10 - \
| samtools sort -@ 8 - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_trimmed.fastq.gz}.sorted.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done")

# Command line entry
### inDir=/Users/coke6162/project20xx_chipseq/trimmomatic \
### outDir=/Users/coke6162/project20xx_chipseq/bwa \
### sbatch --array 0-8 /Users/coke6162/project20xx_chipseq/scripts/bwa.sbatch
