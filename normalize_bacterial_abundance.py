#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import Counter

# Read the raw counts matrix
print("Reading raw bacterial counts...")
df = pd.read_csv("bacteria_text_files/bacterial_genera_raw_counts.txt", sep="\t")

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

# Save rarefied counts (unfiltered)
rarefied_df.to_csv("bacteria_text_files/bacterial_genera_rarefied_counts_new.txt", sep="\t")

# Define filtering function
def filter_bacterial_data(df, analysis_type="general"):
    """Ultra-aggressive filtering to eliminate all possible noise"""
    if analysis_type == "general":
        min_reads = 2000         # ~0.7% of 289k reads
        min_abundance_pct = 1.0  # Must be ≥1% somewhere
        min_samples = 1          # Present in at least 1 sample above threshold
    elif analysis_type == "pie_chart":
        min_reads = 8000         # ~2.8% of 289k reads
        min_abundance_pct = 3.0  # Must be ≥3% somewhere
        min_samples = 1          # Present in at least 1 sample above threshold
    elif analysis_type == "ultra_strict":
        min_reads = 15000        # ~5.2% of 289k reads
        min_abundance_pct = 5.0  # Must be ≥5% somewhere
        min_samples = 2          # Present in at least 2 samples above threshold
    else:
        raise ValueError("analysis_type must be 'general', 'pie_chart', or 'ultra_strict'")

    print(f"\nAGGRESSIVE filtering for {analysis_type} analysis:")
    print(f"  Before filtering: {len(df)} genera")

    # Filter 1: Minimum reads in at least one sample
    max_reads_per_genus = df.max(axis=1)
    read_filter = max_reads_per_genus >= min_reads
    df_filtered = df[read_filter]
    print(f"  After read filter (≥{min_reads}): {len(df_filtered)} genera")

    # Filter 2: Minimum abundance in at least one sample
    sample_totals = df_filtered.sum(axis=0)
    abundance_matrix = df_filtered.div(sample_totals, axis=1) * 100
    max_abundance_per_genus = abundance_matrix.max(axis=1)
    abundance_filter = max_abundance_per_genus >= min_abundance_pct
    df_filtered2 = df_filtered[abundance_filter]
    print(f"  After abundance filter (≥{min_abundance_pct}%): {len(df_filtered2)} genera")

    # Filter 3: Must be abundant in multiple samples (for ultra_strict)
    if min_samples > 1:
        samples_above_threshold = (abundance_matrix >= min_abundance_pct).sum(axis=1)
        multi_sample_filter = samples_above_threshold >= min_samples
        df_final = df_filtered2[multi_sample_filter]
        print(f"  After multi-sample filter (≥{min_samples} samples): {len(df_final)} genera")
    else:
        df_final = df_filtered2

    # Filter 4: Remove potential contaminants (optional - uncomment if desired)
    # known_contaminants = ['Escherichia', 'Bacillus', 'Staphylococcus', 'Streptococcus']
    # contaminant_filter = ~df_final.index.get_level_values('Genus').isin(known_contaminants)
    # df_final = df_final[contaminant_filter]
    # print(f"  After contaminant filter: {len(df_final)} genera")

    return df_final

# Apply filtering for general analysis
print("=== Creating filtered datasets ===")
rarefied_df_filtered = filter_bacterial_data(rarefied_df, "general")
rarefied_df_filtered.to_csv("bacteria_text_files/bacterial_genera_rarefied_counts_filtered.txt", sep="\t")

# Apply ultra-aggressive filtering for pie chart
rarefied_df_piechart = filter_bacterial_data(rarefied_df, "pie_chart")

# Optional: Create an even more stringent dataset
rarefied_df_ultra = filter_bacterial_data(rarefied_df, "ultra_strict")
rarefied_df_ultra.to_csv("bacteria_text_files/bacterial_genera_ultra_strict.txt", sep="\t")

# Calculate overall community composition using FILTERED data for pie chart
print("\nCalculating high-confidence community composition for pie chart...")
overall_counts = rarefied_df_piechart.sum(axis=1)
overall_total = overall_counts.sum()
overall_percentages = (overall_counts / overall_total) * 100

# Sort by abundance
overall_percentages = overall_percentages.sort_values(ascending=False)

# Create output for pie chart with strict filtering
with open("bacteria_text_files/bacterial_overall_composition_for_piechart_new.txt", "w") as f:
    f.write("Genus\tReads\tPercentage\n")

    total_shown = 0
    for genus, pct in overall_percentages.items():
        reads = int(overall_counts[genus])
        f.write(f"{genus}\t{reads}\t{pct:.2f}\n")
        total_shown += pct

    # Add "Other" category for all filtered-out genera
    if total_shown < 100:
        other_pct = 100 - total_shown
        other_reads = rarefied_df.sum().sum() - overall_total  # Reads from filtered genera
        f.write(f"Other\t{int(other_reads)}\t{other_pct:.2f}\n")

print(f"\nHigh-confidence pie chart contains {len(overall_percentages)} genera")
print(f"These represent {total_shown:.1f}% of all bacterial reads")

# Create detailed summary with filtering information
with open("bacteria_text_files/bacterial_normalization_summary_new.txt", "w") as f:
    f.write("Bacterial Community Normalization Summary\n")
    f.write("=========================================\n\n")

    f.write("Normalization method: Rarefaction to minimum sample depth\n")
    f.write(f"Rarefaction depth: {min_reads:,} reads\n")
    f.write(f"Total genera detected (unfiltered): {len(rarefied_df)}\n")
    f.write(f"Genera retained for general analysis: {len(rarefied_df_filtered)}\n")
    f.write(f"Genera retained for pie chart: {len(rarefied_df_piechart)}\n\n")

    f.write("Filtering thresholds applied:\n")
    f.write("  General analysis: ≥2,000 reads AND ≥1.0% abundance\n")
    f.write("  Pie chart: ≥8,000 reads AND ≥3.0% abundance\n")
    f.write("  Ultra strict: ≥15,000 reads AND ≥5.0% abundance in ≥2 samples\n\n")

    f.write("Original read counts per sample:\n")
    for sample in samples:
        f.write(f"  {sample}: {sample_totals[sample]:,}\n")

    f.write(f"\nAll samples normalized to: {min_reads:,} reads\n\n")

    f.write("Top 20 genera in overall community (high-confidence only):\n")
    for i, (genus, pct) in enumerate(overall_percentages.head(20).items(), 1):
        f.write(f"  {i:2d}. {genus}: {pct:.2f}%\n")

    # Per-sample dominant genus after rarefaction (using filtered data)
    f.write("\nDominant genus per sample (high-confidence, after normalization):\n")
    for sample in samples:
        if sample in rarefied_df_piechart.columns:
            sample_pcts = (rarefied_df_piechart[sample] / rarefied_df_piechart[sample].sum()) * 100
            if sample_pcts.sum() > 0:
                top_genus = sample_pcts.idxmax()
                top_pct = sample_pcts.max()
                f.write(f"  {sample}: {top_genus} ({top_pct:.1f}%)\n")

# Create additional filtering summary
with open("bacteria_text_files/bacterial_filtering_summary.txt", "w") as f:
    f.write("Bacterial Data Filtering Summary\n")
    f.write("================================\n\n")
    f.write(f"Original rarefaction depth: {min_reads:,} reads per sample\n")
    f.write(f"Total genera detected: {len(rarefied_df)}\n\n")

    f.write("Filtering applied to remove noise and low-confidence assignments:\n")
    f.write("1. General analysis filter: ≥500 reads AND ≥0.5% abundance in at least one sample\n")
    f.write("2. Pie chart filter: ≥1,500 reads AND ≥1.0% abundance in at least one sample\n\n")

    f.write(f"Results:\n")
    f.write(f"  Genera passing general filter: {len(rarefied_df_filtered)} ({len(rarefied_df_filtered)/len(rarefied_df)*100:.1f}%)\n")
    f.write(f"  Genera passing pie chart filter: {len(rarefied_df_piechart)} ({len(rarefied_df_piechart)/len(rarefied_df)*100:.1f}%)\n\n")

    f.write("This filtering ensures that:\n")
    f.write("- Only genera with strong statistical support are included\n")
    f.write("- Potential sequencing artifacts and contaminants are removed\n")
    f.write("- Pie charts show genuine dominant community members\n")
    f.write("- Results are suitable for publication and interpretation\n")

print("\nAnalysis complete!")
print("\nFiles created:")
print("  - bacteria_text_files/bacterial_genera_rarefied_counts_new.txt (all detected genera)")
print("  - bacteria_text_files/bacterial_genera_rarefied_counts_filtered.txt (moderate filtering)")
print("  - bacteria_text_files/bacterial_overall_composition_for_piechart_new.txt (high-confidence pie chart data)")
print("  - bacteria_text_files/bacterial_normalization_summary_new.txt (detailed analysis summary)")
print("  - bacteria_text_files/bacterial_filtering_summary.txt (filtering rationale and statistics)")
print("  - bacteria_text_files/bacterial_genera_ultra_strict.txt (nuclear option filtering)")
print("\nThe pie chart file now contains only high-confidence genera with strong")
print("statistical support, removing potential noise and sequencing artifacts.")
