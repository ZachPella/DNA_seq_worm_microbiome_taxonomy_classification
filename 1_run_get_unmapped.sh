#!/bin/bash
#SBATCH --job-name=unmapped_for_assembly
#SBATCH --output=unmapped_assembly_%A_%a.out
#SBATCH --error=unmapped_assembly_%A_%a.err
#SBATCH --array=1-9
#SBATCH --mem=15G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1

# Change to the working directory
cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

module load samtools

# Create an array of the BAM files
BAM_FILES=($(ls *.sorted.refrename.bam))

# Get the current BAM file based on the array task ID
INDEX=$((SLURM_ARRAY_TASK_ID-1))
CURRENT_BAM=${BAM_FILES[$INDEX]}

# Extract the sample name without the extension
SAMPLE_NAME=$(basename "${CURRENT_BAM}" .sorted.refrename.bam)

echo "Processing ${SAMPLE_NAME}..."

# Create output directory if it doesn't exist
mkdir -p fastq_unmapped_for_taxonomy_id

# Extract read pairs where both ends are unmapped and convert directly to FASTQ
echo "Extracting unmapped paired reads to FASTQ..."
samtools view -f 12 -F 256 -u "${CURRENT_BAM}" | \
  samtools sort -n | \
  samtools fastq -1 fastq_for_assembly/"${SAMPLE_NAME}.unmapped.R1.fastq" \
                 -2 fastq_for_assembly/"${SAMPLE_NAME}.unmapped.R2.fastq" \
                 -s fastq_for_assembly/"${SAMPLE_NAME}.unmapped.singleton.fastq" \
                 -

echo "Completed processing ${SAMPLE_NAME}"
