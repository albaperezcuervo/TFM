# =============================
# Required Libraries
# =============================
# Load necessary libraries for the batch correction analysis
library(dada2)
library(MMUPHin)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ape)

# Extract metadata from phyloseq objects and merge them
merged_metadata <- rbind(
  data.frame(physeq_16s@sam_data) %>% select(diagnosis, Study),
  data.frame(physeq_wgs@sam_data) %>% select(diagnosis, Study)
)

# Merge OTU tables and transpose
merged_otu_table <- t(dada2::mergeSequenceTables(
  t(physeq_16s@otu_table), 
  t(physeq_wgs@otu_table)
))

# Apply batch effect correction
fit_adjust_batch <- adjust_batch(
  feature_abd = merged_otu_table,
  batch = "Study",
  data = merged_metadata,
  control = list(verbose = TRUE)
)

# Extract metadata from phyloseq object
merged_metadata <- data.frame(sample_data(phyl_div))
merged_metadata <- merged_metadata %>% select(diagnosis, Study, Sex)

# Transpose OTU table to match MMUPHin format
merged_otu_table <- data.frame(t(otu_table(phyl_div)))
corrected_otu_table <- t(merged_otu_table)

# Reapply MMUPHin batch effect correction
fit_adjust_batch <- adjust_batch(
  feature_abd = corrected_otu_table,
  batch = "Study",  # Use 'Study' as batch variable
  data = merged_metadata,
  control = list(verbose = TRUE)
)

# Generate a new phyloseq object
adjusted_otu_table <- fit_adjust_batch$feature_abd_adj
source("https://raw.githubusercontent.com/waldronlab/presentations/master/Waldron_2016-06-07_EPIC/metaphlanToPhyloseq.R")

adjusted_phyloseq <- metaphlanToPhyloseq(
  adjusted_otu_table,
  metadat = merged_metadata,
  simplenames = FALSE
)

# Generate a random phylogenetic tree
random_tree <- rtree(
  ntaxa(adjusted_phyloseq), rooted = TRUE, 
  tip.label = taxa_names(adjusted_phyloseq)
)

# Merge tree with phyloseq object
adjusted_phyloseq <- merge_phyloseq(adjusted_phyloseq, random_tree)

# Save the corrected phyloseq object
saveRDS(adjusted_phyloseq, file = "adjusted_phyloseq.rds")

# Custom function for enhancing plots
custom_plot <- function(plot_obj, title_text) {
  plot_obj +
    geom_point(size = 3, alpha = 0.7, stroke = 0.8) +  
    scale_color_manual(
      name = "Dataset",  # Change legend title
      values = c("0" = "#FDB462", "1" = "#5A9BD4"),  
      labels = c("0" = "16S", "1" = "WGS")  
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 40, face = "bold"), 
      axis.title.y = element_text(size = 40, face = "bold"), 
      axis.text = element_text(size = 35), 
      plot.title = element_text(size = 40, face = "bold", hjust = 0.5),  
      legend.text = element_text(size = 30),  
      legend.title = element_text(size = 35, face = "bold"),  
      legend.key.size = unit(2, "cm")  
    ) +
    ggtitle(title_text) 
}

# Apply compositional transformation to the phyloseq object
pseq.compositional <- transform(adjusted_phyloseq, "compositional")

# Compute beta-diversity using Unifrac distance
ordinated_taxa_unifrac <- ordinate(
  pseq.compositional, method="MDS", distance="bray", weighted=TRUE
)

# Plot ordination for different axis combinations
p1 <- plot_ordination(pseq.compositional, ordinated_taxa_unifrac, color = "Study", axes = c(1, 2)) 
p1 <- custom_plot(p1, "Unifrac MDS Analysis (Axis 1 and 2)")
p1

p2 <- plot_ordination(pseq.compositional, ordinated_taxa_unifrac, color = "Study", axes = c(1, 3)) 
p2 <- custom_plot(p2, "Unifrac MDS Analysis (Axis 1 and 3)")
p2

p3 <- plot_ordination(pseq.compositional, ordinated_taxa_unifrac, color = "Study", axes = c(2, 3)) 
p3 <- custom_plot(p3, "Unifrac MDS Analysis (Axis 2 and 3)")
p3

# Save final corrected phyloseq object
saveRDS(phyl_div, file = "phyl_div_corrected.rds")