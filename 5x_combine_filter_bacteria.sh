#!/bin/bash
#SBATCH --job-name=combine_filter_bacteria
#SBATCH --output=combine_filter_bacteria_%j.out
#SBATCH --error=combine_filter_bacteria_%j.err
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2

cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

echo "Starting bacterial genera extraction and combination..."

# Process each sample
for kraken_report in kraken2_output/*.kraken2.report; do
    if [ -f "$kraken_report" ]; then
        sample=$(basename "$kraken_report" .kraken2.report)
        echo "Processing $sample..."

        # Extract bacterial genus taxids from this sample's Kraken report
        awk '/D.*2.*Bacteria/{flag=1} /D.*2759.*Eukaryota|D.*2157.*Archaea|D.*10239.*Viruses/{flag=0} flag' "$kraken_report" | \
        grep -E "G[[:space:]]+" | awk '{print $5}' > "${sample}_bacterial_taxids.tmp"

        # Filter the Bracken genus output for just bacterial genera
        bracken_file="bracken_output/${sample}.bracken.G.output"
        if [ -f "$bracken_file" ]; then
            # Header
            head -1 "$bracken_file" > "bracken_output/${sample}.bacteria_only.tmp"

            # Filter for bacterial entries
            while read taxid; do
                grep -P "^\S+\t${taxid}\t" "$bracken_file"
            done < "${sample}_bacterial_taxids.tmp" >> "bracken_output/${sample}.bacteria_only.tmp"

            echo "  Found $(( $(wc -l < "bracken_output/${sample}.bacteria_only.tmp") - 1 )) bacterial genera"
        fi

        # Clean up
        rm -f "${sample}_bacterial_taxids.tmp"
    fi
done

# Now combine all the filtered files
echo -e "\nCombining bacterial data..."

cat > combine_bacterial_data.py << 'EOF'
#!/usr/bin/env python3
import os
import glob
from collections import defaultdict

# Get all bacterial filtered files
bacteria_files = sorted(glob.glob("bracken_output/*.bacteria_only.tmp"))

if not bacteria_files:
    print("ERROR: No bacterial filtered files found")
    exit(1)

print(f"Found {len(bacteria_files)} samples with bacterial data")

# Store all data
genus_data = defaultdict(lambda: defaultdict(int))
sample_names = []
genus_to_taxid = {}

# Process each file
for bf in bacteria_files:
    # Extract sample name
    sample = os.path.basename(bf).replace(".bacteria_only.tmp", "")
    sample_names.append(sample)

    with open(bf, 'r') as f:
        header = f.readline()  # Skip header
        count = 0

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                genus_name = parts[0]
                taxid = parts[1]
                reads = int(parts[5])  # new_est_reads

                genus_to_taxid[genus_name] = taxid
                genus_data[genus_name][sample] = reads
                count += 1

        print(f"  {sample}: {count} bacterial genera")

# Write RAW COUNTS matrix
print("\nWriting output files...")
with open("bacterial_genera_raw_counts.txt", 'w') as out:
    out.write("Genus\tTaxID\t" + "\t".join(sample_names) + "\n")

    for genus in sorted(genus_data.keys()):
        out.write(f"{genus}\t{genus_to_taxid[genus]}")
        for sample in sample_names:
            out.write(f"\t{genus_data[genus][sample]}")
        out.write("\n")

# Calculate totals for relative abundance
sample_totals = {s: sum(genus_data[g][s] for g in genus_data) for s in sample_names}

# Write RELATIVE ABUNDANCE matrix (percentages)
with open("bacterial_genera_relative_abundance.txt", 'w') as out:
    out.write("Genus\tTaxID\t" + "\t".join(sample_names) + "\n")

    for genus in sorted(genus_data.keys()):
        out.write(f"{genus}\t{genus_to_taxid[genus]}")
        for sample in sample_names:
            if sample_totals[sample] > 0:
                pct = (genus_data[genus][sample] / sample_totals[sample]) * 100
                out.write(f"\t{pct:.4f}")
            else:
                out.write("\t0.0000")
        out.write("\n")

# Summary
print(f"\nTotal bacterial genera across all samples: {len(genus_data)}")
print("Total bacterial reads per sample:")
for sample in sample_names:
    print(f"  {sample}: {sample_totals[sample]:,}")

# Write summary
with open("bacterial_matrix_summary.txt", 'w') as out:
    out.write("Bacterial Community Analysis Summary\n")
    out.write("====================================\n\n")
    out.write(f"Total bacterial genera: {len(genus_data)}\n")
    out.write(f"Number of samples: {len(sample_names)}\n\n")

    out.write("Bacterial reads per sample:\n")
    for sample in sample_names:
        out.write(f"  {sample}: {sample_totals[sample]:,}\n")

    # Top genera by average abundance
    out.write("\nTop 15 genera by average relative abundance:\n")
    genus_avg = {}
    for genus in genus_data:
        total_pct = sum((genus_data[genus][s] / sample_totals[s]) * 100
                       for s in sample_names if sample_totals[s] > 0)
        genus_avg[genus] = total_pct / len([s for s in sample_names if sample_totals[s] > 0])

    for i, (genus, avg_pct) in enumerate(sorted(genus_avg.items(),
                                               key=lambda x: x[1], reverse=True)[:15], 1):
        out.write(f"  {i:2d}. {genus}: {avg_pct:.2f}%\n")

print("\nOutput files created successfully!")
EOF

python3 combine_bacterial_data.py

# Clean up temporary files
rm -f bracken_output/*.bacteria_only.tmp

echo -e "\nDone! Files created:"
echo "  - bacterial_genera_raw_counts.txt"
echo "  - bacterial_genera_relative_abundance.txt"
echo "  - bacterial_matrix_summary.txt"
