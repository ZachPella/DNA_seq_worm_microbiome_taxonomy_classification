#!/bin/bash
#SBATCH --job-name=normalize_bacteria
#SBATCH --output=normalize_bacteria_%j.out
#SBATCH --error=normalize_bacteria_%j.err
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1

cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

echo "Creating normalized bacterial abundance analysis..."

cat > normalize_bacterial_abundance.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import Counter

# Read the raw counts matrix
print("Reading raw bacterial counts...")
df = pd.read_csv("bacterial_genera_raw_counts.txt", sep="\t")

# Set genus and taxid as index
df = df.set_index(['Genus', 'TaxID'])

# Get sample names
samples = df.columns.tolist()
print(f"Found {len(samples)} samples")

# Calculate total reads per sample
sample_totals = df.sum(axis=0)
print("\nRaw bacterial reads per sample:")
for sample in samples:
    print(f"  {sample}: {sample_totals[sample]:,}")

# Find minimum read count for rarefaction
min_reads = int(sample_totals.min())
print(f"\nMinimum read count: {min_reads:,}")
print("Will rarefy all samples to this depth for fair comparison")

# Perform rarefaction (random subsampling to minimum depth)
print("\nPerforming rarefaction...")
rarefied_data = {}

np.random.seed(42)  # For reproducibility

for sample in samples:
    print(f"  Rarefying {sample}...")

    # Get all reads for this sample
    genus_counts = df[sample]
    genus_counts = genus_counts[genus_counts > 0]  # Remove zeros

    # Create a pool of all reads (just genus names, not the full index)
    read_pool = []
    for (genus, taxid), count in genus_counts.items():
        read_pool.extend([genus] * int(count))

    # Randomly sample to minimum depth
    sampled_reads = np.random.choice(read_pool, size=min_reads, replace=False)

    # Count occurrences
    from collections import Counter
    rarefied_counts = Counter(sampled_reads)
    rarefied_data[sample] = rarefied_counts

# Create rarefied dataframe
print("\nCreating rarefied matrix...")
all_genera = set()
for sample_data in rarefied_data.values():
    all_genera.update(sample_data.keys())

rarefied_df = pd.DataFrame(index=sorted(all_genera), columns=samples)
rarefied_df = rarefied_df.fillna(0)

for sample, counts in rarefied_data.items():
    for genus, count in counts.items():
        rarefied_df.loc[genus, sample] = count

# Save rarefied counts
rarefied_df.to_csv("bacterial_genera_rarefied_counts.txt", sep="\t")

# Calculate overall community composition (sum across all samples)
print("\nCalculating overall community composition...")
overall_counts = rarefied_df.sum(axis=1)
overall_total = overall_counts.sum()
overall_percentages = (overall_counts / overall_total) * 100

# Sort by abundance
overall_percentages = overall_percentages.sort_values(ascending=False)

# Create output for pie chart
with open("bacterial_overall_composition_for_piechart.txt", "w") as f:
    f.write("Genus\tReads\tPercentage\n")

    # Write top genera individually
    cumulative_pct = 0
    other_count = 0
    other_pct = 0

    for genus, pct in overall_percentages.items():
        reads = int(overall_counts[genus])
        if cumulative_pct < 90 and pct >= 0.5:  # Top genera until 90% or >= 0.5%
            f.write(f"{genus}\t{reads}\t{pct:.2f}\n")
            cumulative_pct += pct
        else:
            other_count += reads
            other_pct += pct

    # Add "Other" category
    if other_count > 0:
        f.write(f"Other\t{int(other_count)}\t{other_pct:.2f}\n")

print(f"\nWrote pie chart data with {len(overall_percentages)} genera")
if cumulative_pct > 0:
    print(f"Top genera account for {cumulative_pct:.1f}% of reads")
if other_pct > 0:
    print(f"'Other' category contains {other_pct:.1f}% of reads")

# Also create a detailed summary
with open("bacterial_normalization_summary.txt", "w") as f:
    f.write("Bacterial Community Normalization Summary\n")
    f.write("=========================================\n\n")

    f.write("Normalization method: Rarefaction to minimum sample depth\n")
    f.write(f"Rarefaction depth: {min_reads:,} reads\n")
    f.write(f"Total genera detected: {len(overall_percentages)}\n\n")

    f.write("Original read counts per sample:\n")
    for sample in samples:
        f.write(f"  {sample}: {sample_totals[sample]:,}\n")

    f.write(f"\nAll samples normalized to: {min_reads:,} reads\n\n")

    f.write("Top 20 genera in overall community:\n")
    for i, (genus, pct) in enumerate(overall_percentages.head(20).items(), 1):
        f.write(f"  {i:2d}. {genus}: {pct:.2f}%\n")

    # Per-sample dominant genus after rarefaction
    f.write("\nDominant genus per sample (after normalization):\n")
    for sample in samples:
        sample_pcts = (rarefied_df[sample] / rarefied_df[sample].sum()) * 100
        top_genus = sample_pcts.idxmax()
        top_pct = sample_pcts.max()
        f.write(f"  {sample}: {top_genus} ({top_pct:.1f}%)\n")

print("\nAnalysis complete!")
print("\nFiles created:")
print("  - bacterial_genera_rarefied_counts.txt (normalized count matrix)")
print("  - bacterial_overall_composition_for_piechart.txt (ready for pie chart)")
print("  - bacterial_normalization_summary.txt (analysis summary)")
print("\nThe pie chart file contains the overall bacterial composition across all samples,")
print("with proper normalization for different sequencing depths.")
EOF

# Run the normalization script
python3 normalize_bacterial_abundance.py

# Show the pie chart data
echo -e "\n=== Pie chart data (top entries) ==="
head -20 bacterial_overall_composition_for_piechart.txt

echo -e "\n=== Summary ==="
tail -20 bacterial_normalization_summary.txt
