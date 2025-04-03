# =============================
# Required Libraries
# =============================
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


# =============================
# Assess Normality 
# =============================


abund_table <- as.data.frame(otu_table(pseq.filtered))
metadata <- as.data.frame(sample_data(pseq.filtered))


shapiro_results <- list()


for (taxa in rownames(abund_table)) {
  taxa_abundance <- abund_table[taxa, ]
  df <- data.frame(
    abundance = as.numeric(taxa_abundance),
    diagnosis = metadata$diagnosis
  )
  control <- df$abundance[df$diagnosis == "Control"]
  alzheimer <- df$abundance[df$diagnosis == "Alzheimer"]
  shapiro_results[[taxa]] <- list(
    control_p_value = shapiro.test(control)$p.value,
    alzheimer_p_value = shapiro.test(alzheimer)$p.value
  )
}


shapiro_df <- do.call(rbind, lapply(names(shapiro_results), function(taxa) {
  data.frame(
    taxa = taxa,
    control_p_value = shapiro_results[[taxa]]$control_p_value,
    alzheimer_p_value = shapiro_results[[taxa]]$alzheimer_p_value
  )
}))


shapiro_df$control_normal <- ifelse(shapiro_df$control_p_value > 0.05, "Normal", "Non-Normal")
shapiro_df$alzheimer_normal <- ifelse(shapiro_df$alzheimer_p_value > 0.05, "Normal", "Non-Normal")
print(shapiro_df)


# =============================
# Relative Abundance by Diagnosis - PHYLUM
# =============================


tb_phylum <- psmelt(pseq.filtered)


# Get Top 10 Most Abundant Phyla
top10_phyla <- tb_phylum %>%
  group_by(Phylum) %>%
  summarise(median_abund = median(Abundance)) %>%
  arrange(desc(median_abund)) %>%
  slice_head(n = 10)


tbphylum <- tb_phylum %>% filter(Phylum %in% top10_phyla$Phylum)
tbphylum$diagnosis <- factor(tbphylum$diagnosis, levels = c("Control", "Alzheimer"))


# Barplot
plot_phylum <- ggplot(tbphylum, aes(diagnosis, Abundance, fill = Phylum)) + 
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  theme_light() +
  xlab("Diagnosis") + ylab("Relative Abundance") +
  ggtitle("Relative Abundance of Top 10 Phyla by Diagnosis") +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 25, face = "bold")
  )
print(plot_phylum)


# Wilcoxon Test
wilcoxon_results <- tbphylum %>%
  group_by(Phylum) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ diagnosis, data = cur_data())$p.value
  ) %>%
  mutate(significant = p_value < 0.05)


print(wilcoxon_results)


# Filter Significant Phyla
tbphylum_significant <- tbphylum %>%
  filter(Phylum %in% wilcoxon_results$Phylum[wilcoxon_results$significant == TRUE])


# Violin Plot for Significant Phyla
violinplot_phyla <- ggplot(tbphylum_significant, aes(x = diagnosis, y = log10(Abundance + 1e-6), fill = diagnosis)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, color = "black") +
  geom_jitter(aes(color = diagnosis), position = position_jitter(width = 0.2), alpha = 0.5, size = 1, show.legend = FALSE) +
  facet_wrap(~Phylum, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 7, fontface = "bold", label.y.npc = "top") +
  scale_fill_manual(values = c("Control" = "#B3DE69", "Alzheimer" = "#FA8072")) +
  scale_color_manual(values = c("Control" = "#B3DE69", "Alzheimer" = "#FA8072")) +
  theme_minimal() +
  labs(
    title = "Relative Abundance of Statistically Significant Phyla",
    x = "Diagnosis",
    y = "Log10 Relative Abundance"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )
print(violinplot_phyla)


# By Sex (Facet)
plot_phylum_sex <- ggplot(tbphylum, aes(diagnosis, Abundance, fill = Phylum)) +
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Sex) +
  theme_light() +
  xlab("Diagnosis") + ylab("Relative Abundance") +
  ggtitle("Relative Abundance of Top 10 Phyla by Diagnosis and Sex") +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1.5, "cm"),
    strip.text = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  )
print(plot_phylum_sex)


# ANOVA
anova_results <- aov(Abundance ~ Sex * diagnosis, data = tbphylum)
summary(anova_results)


# =============================
# Repeat for GENUS Level
# =============================


tb_genus <- psmelt(pseq.filtered)


# Top 10 Most Abundant Genera
top10_genus <- tb_genus %>%
  group_by(Genus) %>%
  summarise(median_abund = median(Abundance)) %>%
  arrange(desc(median_abund)) %>%
  slice_head(n = 10)


tbgenus <- tb_genus %>% filter(Genus %in% top10_genus$Genus)
tbgenus$diagnosis <- factor(tbgenus$diagnosis, levels = c("Control", "Alzheimer"))


# Barplot
plot_genus <- ggplot(tbgenus, aes(diagnosis, Abundance, fill = Genus)) + 
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Paired") +
  theme_light() +
  xlab("Diagnosis") + ylab("Relative Abundance") +
  ggtitle("Relative Abundance of Top 10 Genera by Diagnosis") +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 25, face = "bold")
  )
print(plot_genus)


# Wilcoxon Test for Genera
wilcoxon_genus_results <- tbgenus %>%
  group_by(Genus) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ diagnosis, data = cur_data())$p.value
  ) %>%
  mutate(significant = p_value < 0.05)


print(wilcoxon_genus_results)


# Filter Significant Genera
tbgenus_significant <- tbgenus %>%
  filter(Genus %in% wilcoxon_genus_results$Genus[wilcoxon_genus_results$significant == TRUE])


# Violin Plot for Significant Genera
violinplot_genus <- ggplot(tbgenus_significant, aes(x = diagnosis, y = log10(Abundance + 1e-6), fill = diagnosis)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, color = "black") +
  geom_jitter(aes(color = diagnosis), position = position_jitter(width = 0.2), alpha = 0.5, size = 1, show.legend = FALSE) +
  facet_wrap(~Genus, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 7, fontface = "bold", label.y.npc = "top") +
  scale_fill_manual(values = c("Control" = "#B3DE69", "Alzheimer" = "#FA8072")) +
  scale_color_manual(values = c("Control" = "#B3DE69", "Alzheimer" = "#FA8072")) +
  theme_minimal() +
  labs(
    title = "Relative Abundance of Statistically Significant Genera",
    x = "Diagnosis",
    y = "Log10 Relative Abundance"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )
print(violinplot_genus)


# By Sex (Facet)

plot_genus_sex <- ggplot(tbgenus, aes(diagnosis, Abundance, fill = Genus)) + 
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Paired") +
  facet_wrap(~Sex) +
  theme_light() +
  xlab("Diagnosis") + ylab("Relative Abundance") +
  ggtitle("Relative Abundance of Top 10 Genera by Diagnosis and Sex") +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1.5, "cm"),
    strip.text = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  )
print(plot_genus_sex)

# ANOVA

anova_genus_results <- aov(Abundance ~ Sex * diagnosis, data = tbgenus)
summary(anova_genus_results)


