# =============================
# Required Libraries
# =============================
# Load necessary libraries for microbiome analysis, visualization, and data manipulation
library(tidyverse)
library(ape)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(MMUPHin)
library(ggpubr)
library(vegan)
library(genefilter)
library(microViz)
library(ggtree)
library(tibble)
library(ggridges)
library(taxizedb)
library(data.table)

# Load phyloseq object at genus level
physeq.all.g <- readRDS("physeq_all_tree_Genus.rds")

# Calculate and visualize sequencing depth per sample
readcount <- data.table(as(sample_data(physeq.all.g), "data.frame"),
                        TotalReads = sample_sums(physeq.all.g), 
                        keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID2")
ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

# Prepare OTU table for rarefaction curves
otu.rare <- otu_table(physeq.all.g)
otu.rare <- as.data.frame(t(otu.rare))
sample_names <- rownames(otu.rare)

# Generate rarefaction curves to assess sequencing depth adequacy
library(vegan)
otu.rarecurve <- rarecurve(otu.rare, step = 10000, label = TRUE)

# Define colors ensuring Study 0 is the first group
study_group <- sample_data(physeq.all.g)$Study
colors_study <- ifelse(study_group == "0", "#FDB462", "#5A9BD4")  # Orange for Study 0, Blue for Study 1

# Plot rarefaction curves for Study 1 (Blue)
rarecurve(otu.rare[study_group == "1", ], step = 10000, label = FALSE, col = "#5A9BD4", 
          lwd = 2, cex.axis = 1.4, cex.lab = 1.6)

# Plot rarefaction curves for Study 0 (Orange)
rarecurve(otu.rare[study_group == "0", ], step = 10000, label = FALSE, col = "#FDB462", 
          lwd = 2, cex.axis = 1.4, cex.lab = 1.6)

# Determine minimum read count for rarefaction
min_reads <- min(readcount$TotalReads)
min_reads

# Perform rarefaction to equalize sequencing depth across samples
phyl_div <- rarefy_even_depth(physeq.all.g, sample.size = min_reads, rngseed = 42)

# Generate rarefaction curves after rarefaction
otu.rare <- otu_table(phyl_div)
otu.rare <- as.data.frame(t(otu.rare))
sample_names <- rownames(otu.rare)
otu.rarecurve <- rarecurve(otu.rare, step = 1000, label = TRUE)

# Save post-rarefaction plot with study-based color differentiation
study_group <- sample_data(phyl_div)$Study
colors_study <- ifelse(study_group == "0", "#FDB462", "#5A9BD4")
rarecurve(otu.rare, step = 1000, label = FALSE, col = colors_study, lwd = 2, 
          cex.axis = 1.4, cex.lab = 1.6)

# Save the final rarefied phyloseq object
saveRDS(phyl_div, file = "phyl_rarefied.rds")