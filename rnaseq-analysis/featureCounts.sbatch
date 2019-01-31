#!/bin/bash
#SBATCH -p short # Partition or queue
#SBATCH --job-name=project20xx_rnaseq_featureCounts # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=coke6162@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --cpus-per-task=8         # cpus
#SBATCH --mem=64gb                # Memory limit
#SBATCH --time=1:00:00            # Time limit hrs:min:sec
#SBATCH --output=/Users/coke6162/sbatch_script/output/project20xx_rnaseq_featureCounts_%j.out # Standard output
#SBATCH --error=/Users/coke6162/sbatch_script/error/project20xx_rnaseq_featureCounts_%j.err # Standard error log

# Set constant variables
inDir=/Users/coke6162/project20xx_rnaseq/hisat2
outDir=/Users/coke6162/project20xx_rnaseq/featureCounts
annotationFile=/scratch/Shares/chuong/genomes/mm10/gencode.vM18.annotation.gtf # Change this to suit your needs

# Load subread
module load subread

# Run featureCounts, outputting a single table containing count information for all samples
# Note: this does not include any normalization (as required by DESeq2)
pwd; hostname; date

echo "Starting featureCounts gene-level summarization"

for i in {0..2} # "-s 0" for unstranded, "-s 1" for stranded, "-s 2" for reversely stranded
do
featureCounts -O -T 1 -s $i -t exon -g gene_id -a ${annotationFile} -o ${outDir}/featureCounts_s$i.txt ${inDir}/*.sorted.uniq.bam &
done
wait

echo $(date +"[%b %d %H:%M:%S] Done")
