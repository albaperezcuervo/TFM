# =============================
# Required Libraries
# =============================
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(vegan)


# =============================
# Alpha Diversity Analysis
# =============================


# Compute alpha diversity indices
df_richness <- estimate_richness(phyl_div, measures = c("Chao1", "Simpson", "Shannon"))
df_richness$diagnosis <- sample_data(phyl_div)$diagnosis


# Evaluate normality with Shapiro-Wilk test
shapiro_results <- df_richness %>%
  pivot_longer(cols = c("Chao1", "Simpson", "Shannon"), names_to = "index", values_to = "value") %>%
  group_by(index, diagnosis) %>%
  summarise(
    p_value_shapiro = shapiro.test(value)$p.value,
    normality = ifelse(p_value_shapiro > 0.05, "Normal", "Non-Normal"),
    .groups = "drop"
  )


# Perform statistical test (t-test or Wilcoxon) based on normality
stat_results <- df_richness %>%
  pivot_longer(cols = c("Chao1", "Simpson", "Shannon"), names_to = "index", values_to = "value") %>%
  group_by(index) %>%
  summarise(
    test_applied = ifelse(
      all(shapiro_results$normality[shapiro_results$index == unique(index)] == "Normal"),
      "t-test",
      "wilcox-test"
    ),
    p_value = ifelse(
      test_applied == "t-test",
      t.test(value ~ diagnosis, data = cur_data())$p.value,
      wilcox.test(value ~ diagnosis, data = cur_data())$p.value
    ),
    .groups = "drop"
  )


# Combine results for summary table
combined_results <- shapiro_results %>%
  left_join(stat_results, by = "index") %>%
  select(index, diagnosis, p_value_shapiro, normality, test_applied, p_value)


print(combined_results)


# =============================
# Taxa Filtering
# =============================


# Convert to relative abundance (compositional)
pseq.compositional <- transform(phyl_div, "compositional")


# Number of taxa before filtering
taxa_antes <- ntaxa(pseq.compositional)
cat("Number of taxa before filtering:", taxa_antes, "\n")


# Define filtering thresholds
min_abundance <- 1e-4
min_presence <- 0.20


# Filtering function
filter_taxa <- function(taxon) {
  abundance_filter <- taxon > min_abundance
  presence <- sum(abundance_filter) / length(taxon)
  return(presence >= min_presence)
}


# Apply filtering to taxa
pseq.filtered <- prune_taxa(apply(otu_table(pseq.compositional), 1, filter_taxa), pseq.compositional)


# Number of taxa after filtering
taxa_despues <- ntaxa(pseq.filtered)
cat("Number of taxa after filtering:", taxa_despues, "\n")


# Clean and standardize taxonomic names
tax_table(pseq.filtered) <- apply(tax_table(pseq.filtered), 2, trimws)
tax_table(pseq.filtered)[, "Phylum"] <- gsub("^Firmicutes_.*", "Firmicutes", tax_table(pseq.filtered)[, "Phylum"])
tax_table(pseq.filtered)[, "Phylum"] <- gsub("^Methanobacteriota_A_1229$", "Methanobacteriota", tax_table(pseq.filtered)[, "Phylum"])
tax_table(pseq.filtered)[, "Phylum"] <- gsub("^Desulfobacterota_I$", "Desulfobacterota", tax_table(pseq.filtered)[, "Phylum"])
tax_table(pseq.filtered)[, "Class"] <- gsub("^Clostridia_.*", "Clostridia", tax_table(pseq.filtered)[, "Class"])


# =============================
# Beta Diversity Analysis
# =============================


# Compute distance matrices
bray_dist <- phyloseq::distance(pseq.filtered, method = "bray")
jaccard_dist <- phyloseq::distance(pseq.filtered, method = "jaccard")
unifrac_dist <- phyloseq::distance(pseq.filtered, method = "unifrac")


# Extract metadata
metadata_df <- as.data.frame(sample_data(pseq.filtered))


# Perform PERMANOVA
permanova_bray <- adonis2(bray_dist ~ metadata_df$diagnosis)
permanova_jaccard <- adonis2(jaccard_dist ~ metadata_df$diagnosis)
permanova_unifrac <- adonis2(unifrac_dist ~ metadata_df$diagnosis)


# Display PERMANOVA results
print(permanova_bray)
print(permanova_jaccard)
print(permanova_unifrac)


# Ordination Plots Function

custom_plot <- function(plot_obj, title_text) {
  plot_obj +
    geom_point(size = 3, alpha = 0.7, stroke = 0.8) +
    theme_bw() +
    stat_ellipse(aes(group = diagnosis), linetype = 1, level = 0.95) +
    scale_color_manual(values = c("Control" = "#B3DE69", "Alzheimer" = "#FA8072")) +
    theme(
      axis.title.x = element_text(size = 35, face = "bold"),
      axis.title.y = element_text(size = 35, face = "bold"),
      axis.text = element_text(size = 32),
      plot.title = element_text(size = 33, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 30, face = "bold")
    ) +
    ggtitle(title_text)
}


# Generate Ordination Plots

distances <- list(bray = "bray", jaccard = "jaccard", unifrac = "unifrac")
axes_combinations <- list(c(1, 2), c(1, 3), c(2, 3))


for (dist_name in names(distances)) {
  for (axes in axes_combinations) {
    ordinated <- ordinate(pseq.filtered, method = "MDS", distance = distances[[dist_name]])
    plot <- plot_ordination(pseq.filtered, ordinated, axes = axes, color = "diagnosis")
    plot <- custom_plot(plot, paste0(toupper(dist_name), " MDS Analysis (Axes ", paste(axes, collapse = " and "), ")"))
    print(plot)
  }
}


