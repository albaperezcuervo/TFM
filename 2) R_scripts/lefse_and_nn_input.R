# =============================
# Required Libraries
# =============================
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(janitor)


# =============================
# Section 1: Prepare Input for LEfSe -> MODIFICAR
# =============================


# Extract OTU table and convert to data frame
tax_gen <- as.data.frame(pseq.filtered@otu_table@.Data)
tax_gen <- rownames_to_column(tax_gen, var = "ASV")


# Extract and collapse taxonomy to genus level
tax_info <- as.data.frame(pseq.filtered@tax_table@.Data)
tax_info <- tax_info %>%
  unite(Phylum, Class, Order, Family, Genus, col = genus, sep = "|") %>%
  rownames_to_column(var = "ASV")


# Merge taxonomy with OTU table by ASV
tax_gen <- tax_info %>%
  select(ASV, genus) %>%
  left_join(tax_gen, by = "ASV")


# Remove ASV column (now genus is the ID)
tax_gen$ASV <- NULL


# Transpose the table to have samples as rows
tax_gen <- t(tax_gen)
tax_gen <- tax_gen %>% row_to_names(row_number = 1)
tax_gen <- as.data.frame(tax_gen)
tax_gen <- rownames_to_column(tax_gen, var = "SampleID")


# Prepare Metadata for LEfSe

# Extract sample metadata
metadata_gen <- as.data.frame(pseq.filtered@sam_data)
write_csv(metadata_gen, "metadata_gen.csv")  # Export full metadata for backup


# Read and select only relevant columns for LEfSe
metadata_gen <- read_csv("metadata_gen.csv") %>%
  select(SampleID, diagnosis)


# Merge metadata with genus-level abundance table
tax_gen_f <- left_join(metadata_gen, tax_gen)


# Final LEfSe formatting
lefse <- t(tax_gen_f)
lefse[is.na(lefse)] <- 0
lefse <- rownames_to_column(as.data.frame(lefse))


# Export in LEfSe-compatible format (tab-delimited, no headers)
write_tsv(lefse, "lefse_gen.txt", col_names = FALSE)


# =============================
# Section 2: Prepare Input for Neural Networks
# =============================


# Extract OTU and taxonomy tables
otu <- as.data.frame(pseq.filtered@otu_table@.Data)
tax <- as.data.frame(pseq.filtered@tax_table@.Data)


# Add OTU IDs as columns
otu$OTU_ID <- rownames(otu)
tax$OTU_ID <- rownames(tax)


# Merge OTU abundances and taxonomy
combined <- merge(tax, otu, by = "OTU_ID")


# Remove 'Species' column if it exists
if ("Species" %in% colnames(combined)) {
  combined$Species <- NULL
}


# Export merged OTU + taxonomy table
write.csv(combined, "combined_otu_taxonomy.csv", row.names = FALSE)


