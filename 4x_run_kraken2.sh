#!/bin/bash
#SBATCH --job-name=kraken2_classify
#SBATCH --output=kraken2_%A_%a.out
#SBATCH --error=kraken2_%A_%a.err
#SBATCH --array=1-9
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8

cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

# Load correct modules (use Kraken2 2.0.8-beta for Bracken compatibility)
module load kraken2/2.0.8-beta
module load bracken

# Get list of R1 files (use trimmed if available, otherwise use original)
if [ -d "trimmed_fastq" ]; then
    R1_FILES=($(ls trimmed_fastq/*.trimmed.R1.fastq))
    R1_SUFFIX=".trimmed.R1.fastq"
    R2_SUFFIX=".trimmed.R2.fastq"
else
    R1_FILES=($(ls fastq_for_assembly/*.unmapped.R1.fastq))
    R1_SUFFIX=".unmapped.R1.fastq"
    R2_SUFFIX=".unmapped.R2.fastq"
fi

# Get current file based on array task ID
INDEX=$((SLURM_ARRAY_TASK_ID-1))
R1_FILE=${R1_FILES[$INDEX]}
R2_FILE=${R1_FILE/$R1_SUFFIX/$R2_SUFFIX}

# Extract sample name
SAMPLE=$(basename ${R1_FILE} $R1_SUFFIX)

# Create output directories
mkdir -p kraken2_output
mkdir -p bracken_output

# Run Kraken2
kraken2 --db $KRAKEN2_DB \
    --threads 8 \
    --paired \
    --use-names \
    --classified-out kraken2_output/${SAMPLE}_classified#.fastq \
    --unclassified-out kraken2_output/${SAMPLE}_unclassified#.fastq \
    --memory-mapping \
    --report kraken2_output/${SAMPLE}.kraken2.report \
    --output kraken2_output/${SAMPLE}.kraken2.output \
    ${R1_FILE} ${R2_FILE}

# Run Bracken for abundance estimation at different taxonomic levels
for LEVEL in S G F O C P; do
    bracken -d $KRAKEN2_DB \
        -i kraken2_output/${SAMPLE}.kraken2.report \
        -o bracken_output/${SAMPLE}.bracken.${LEVEL}.output \
        -w bracken_output/${SAMPLE}.bracken.${LEVEL}.report \
        -r 100 \
        -l ${LEVEL} \
        -t 10
done
