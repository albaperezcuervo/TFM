# =============================
# Required Libraries
# =============================
# Load necessary libraries for microbiome analysis and data handling
library(tidyverse)
library(ape)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(MMUPHin)
library(ggpubr)
library(phyloseq)
library(vegan)
library(genefilter)
library(microViz)
library(ggtree)
library(tibble)
library(ggridges)
library(taxizedb)

# Creating phyloseq objects
# First, generate phyloseq objects separately (without phylogenetic tree)

## 1) 16S data
# Load ASV and taxonomy data (Qiime2 outputs)
ASV <- read_qza("./Data/16s.gg2.biom.qza")
taxonomyq <- read_qza("./Data/16s.gg2.taxonomy.qza")

# Load metadata for 16S and rename column for consistency
metadata_16s <- read_csv("./Data/metadata_16s.csv")
metadata_16s <- metadata_16s %>% rename(Sex = sex)

# Format taxonomy table by separating into multiple levels and removing prefixes
taxtableq <- taxonomyq$data %>%
  as_tibble() %>%
  separate(Taxon, sep=";", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
  mutate(across(everything(), ~gsub("[dpcogs]__", "", .)))

# Construct the phyloseq object for 16S
OTUs <- otu_table(ASV$data, taxa_are_rows = TRUE)
TAXq <- tax_table(as.data.frame(taxtableq) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix())
sample_metadata_16s <- sample_data(metadata_16s %>% as.data.frame())
sample_names(sample_metadata_16s) <- metadata_16s$SampleID
physeq_16s <- merge_phyloseq(OTUs, TAXq, sample_metadata_16s)
saveRDS(physeq_16s, file = "physeq_16s.rds")

# Verify counts of Control vs Alzheimer (expecting 30 Alzheimer and 30 Control)
table(physeq_16s@sam_data$diagnosis)

## 2) WGS data
# Load ASV and taxonomy data (Qiime2 outputs)
ASV <- read_qza("./Data/woltka_gg2.biom.qza")
taxonomyq <- read_qza("./Data/woltka_gg2.taxonomy.qza")

# Load metadata for WGS and add diagnosis column
metadata_wgs <- read_csv("./Data/metadata_wgs.csv") %>%
  mutate(diagnosis = case_when(
    amyloid.positive.AF == 0 ~ "Control",
    amyloid.positive.AF == 1 ~ "Alzheimer",
    TRUE ~ NA_character_
  ))

# Format taxonomy table similar to 16S
taxtableq <- taxonomyq$data %>%
  as_tibble() %>%
  separate(Taxon, sep=";", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
  mutate(across(everything(), ~gsub("[dpcogs]__", "", .)))

# Construct the phyloseq object for WGS
OTUs <- otu_table(ASV$data, taxa_are_rows = TRUE)
TAXq <- tax_table(as.data.frame(taxtableq) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix())
sample_metadata_wgs <- sample_data(metadata_wgs %>% as.data.frame())
sample_names(sample_metadata_wgs) <- metadata_wgs$SampleID
physeq_wgs <- merge_phyloseq(OTUs, TAXq, sample_metadata_wgs)
saveRDS(physeq_wgs, file = "physeq_wgs.rds")

# Verify counts of Control vs Alzheimer (expecting 49 Alzheimer and 115 Control)
table(physeq_wgs@sam_data$diagnosis)

## 3) merging phyloseq objects
# Add study identifiers before merging phyloseq objects
sample_data(physeq_16s)$Study <- as.factor("0")
sample_data(physeq_wgs)$Study <- as.factor("1")

# Merge phyloseq objects (both lack phylogenetic trees, contain diagnosis info)
physeq.all <- merge_phyloseq(physeq_16s, physeq_wgs)

# Generate a random phylogenetic tree and merge it with the combined phyloseq object
random_tree <- rtree(ntaxa(physeq.all), rooted = TRUE, tip.label = taxa_names(physeq.all))
physeq.all.tree <- merge_phyloseq(physeq.all, random_tree)
saveRDS(physeq.all.tree, file = "physeq_all_tree.rds")

# Agglomerate data at genus level and save
physeq.all.g <- tax_glom(physeq.all.tree, taxrank = "Genus")
saveRDS(physeq.all.g, file = "physeq_all_tree_Genus.rds")
