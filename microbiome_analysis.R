#!/usr/bin/env Rscript

# Microbiome Statistical Analysis
# Based on methods: NMDS with Bray-Curtis, Kruskal-Wallis with FDR correction

# Load libraries
library(vegan, lib.loc = "~/R_libs")
library(ecodist, lib.loc = "~/R_libs")

# Set working directory
setwd("/work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og")

# Read the combined Bracken results
combined_data <- read.table("combined_bracken_species.txt",
                           header = TRUE,
                           sep = "\t",
                           row.names = 1)

print("Data dimensions:")
print(dim(combined_data))
print("Column names:")
print(colnames(combined_data))

# Select only the numeric count columns (those ending with "_num")
count_columns <- grep("_num$", colnames(combined_data), value = TRUE)
print("Count columns found:")
print(count_columns)

# Extract only the count data
count_data <- combined_data[, count_columns]

# Clean up sample names (remove the suffix)
colnames(count_data) <- gsub("\\.bracken\\.S\\.output_num$", "", colnames(count_data))
colnames(count_data) <- gsub("\\.", "-", colnames(count_data))  # Replace dots with dashes

print("Cleaned data dimensions:")
print(dim(count_data))
print("Sample names:")
print(colnames(count_data))

# Transpose data so samples are rows and species are columns
community_matrix <- t(count_data)

# Convert to numeric matrix (in case there are any character values)
community_matrix <- apply(community_matrix, 2, as.numeric)
rownames(community_matrix) <- colnames(count_data)

# Remove species with zero total abundance
community_matrix <- community_matrix[, colSums(community_matrix, na.rm = TRUE) > 0]

# Remove any rows (samples) with all zeros
community_matrix <- community_matrix[rowSums(community_matrix, na.rm = TRUE) > 0, ]

print(paste("Community matrix dimensions:", nrow(community_matrix), "samples x", ncol(community_matrix), "species"))

# Calculate Bray-Curtis dissimilarity
bray_dist <- vegdist(community_matrix, method = "bray")

# Perform NMDS
set.seed(123)  # for reproducibility
nmds_result <- metaMDS(community_matrix,
                      distance = "bray",
                      k = 2,  # 2 dimensions
                      trymax = 100)

print(paste("NMDS Stress:", round(nmds_result$stress, 3)))
if(nmds_result$stress < 0.2) {
  print("Good fit (stress < 0.2)")
} else {
  print("Warning: Poor fit (stress >= 0.2)")
}

# Create NMDS plot
pdf("NMDS_plot.pdf", width = 8, height = 6)
plot(nmds_result, type = "n", main = "NMDS of Microbial Communities")
points(nmds_result, display = "sites", pch = 19, col = "blue")
text(nmds_result, display = "sites", labels = rownames(community_matrix), pos = 3, cex = 0.8)
dev.off()

# Calculate diversity indices
diversity_indices <- data.frame(
  Sample = rownames(community_matrix),
  Shannon = diversity(community_matrix, index = "shannon"),
  Simpson = diversity(community_matrix, index = "simpson"),
  Richness = specnumber(community_matrix),
  Evenness = diversity(community_matrix, index = "shannon") / log(specnumber(community_matrix))
)

print("Diversity indices:")
print(diversity_indices)

# Write diversity indices to file
write.csv(diversity_indices, "diversity_indices.csv", row.names = FALSE)

# If you have metadata with experimental groups, uncomment and modify:
# metadata <- read.csv("sample_metadata.csv")  # Create this file with sample info
#
# # Kruskal-Wallis test for differences between groups
# kw_shannon <- kruskal.test(Shannon ~ Group, data = diversity_indices)
# print("Kruskal-Wallis test for Shannon diversity:")
# print(kw_shannon)
#
# # Post-hoc Dunn test (requires dunn.test package)
# if(kw_shannon$p.value < 0.05) {
#   library(dunn.test)
#   dunn_result <- dunn.test(diversity_indices$Shannon, diversity_indices$Group, method = "fdr")
#   print("Post-hoc Dunn test results:")
#   print(dunn_result)
# }

# Export community matrix for further analysis
write.csv(community_matrix, "community_matrix.csv")

# Export NMDS coordinates
nmds_coords <- data.frame(
  Sample = rownames(community_matrix),
  NMDS1 = nmds_result$points[,1],
  NMDS2 = nmds_result$points[,2]
)
write.csv(nmds_coords, "nmds_coordinates.csv", row.names = FALSE)

print("Analysis complete! Files created:")
print("- NMDS_plot.pdf")
print("- diversity_indices.csv")
print("- community_matrix.csv")
print("- nmds_coordinates.csv")
