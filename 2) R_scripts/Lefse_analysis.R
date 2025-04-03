# =============================
# Required Libraries
# =============================
# Load necessary libraries for the lefse analysis
library(dplyr)
library(stringr)

# =============================
# Lefse Analysis
# =============================
# Read the .res file
file <- "Data_lefse/lefse_gen.res"  
data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract the first column containing taxonomy
taxonomies <- data[, 1]

# Function to check if genus is present in the taxonomy
extract_genus <- function(taxonomy) {
  levels <- unlist(strsplit(taxonomy, "\\."))  
  if (length(levels) >= 5) { 
    return(levels[5])
  } else {
    return(NA)  
  }
}

# Apply the function and add the "Genus" column to the original dataframe
data$Genus <- sapply(taxonomies, extract_genus)

# Filter rows where Genus is not NA and the fourth column is greater than 2.5
filtered_data <- data %>%
  filter(!is.na(Genus) & data[, 4] > 2.5)



# Display the "Genus" column from the filtered dataset
filtered_data$Genus


# Replace specific names in the Genus column
filtered_data$Genus <- gsub("_B$", "", filtered_data$Genus)  
filtered_data$Genus <- gsub("_A$", "", filtered_data$Genus)  
filtered_data$Genus <- gsub("_[0-9]+$", "", filtered_data$Genus)  

# List of names to be removed
names_to_remove <- c("CAG_1427", "ER4", "CAG")

# Filter out unwanted names in Genus
filtered_data <- filtered_data %>%
  filter(!Genus %in% names_to_remove)

# Display the "Genus" column after cleaning
filtered_data$Genus



# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Rename columns for clarity
colnames(filtered_data)[1] <- "Taxonomy"
colnames(filtered_data)[3] <- "Group"
colnames(filtered_data)[4] <- "LDA_Score"

# Adjust LDA Score values:
# - Make Alzheimer values negative so they appear on the left side of 0
filtered_data <- filtered_data %>%
  mutate(LDA_Score = ifelse(Group == "Alzheimer", -LDA_Score, LDA_Score))

# Define colors for the groups
colors <- c("Alzheimer" = "#FA8072", "Control" = "#B3DE69")



# Create and display the plot
p <- ggplot(filtered_data, aes(x = reorder(Genus, LDA_Score), y = LDA_Score, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Make the plot horizontal
  scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = seq(-4, 4, by = 1), limits = c(-4, 4)) +  
  labs(x = "", y = "LDA SCORE (log 10)", fill = "") +  
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 35, face = "bold"),  
    legend.title = element_text(size = 30, face = "bold"), 
    axis.text.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_text(size = 25, face = "bold"),  
    axis.title.y = element_text(size = 25, face = "bold"), 
    axis.title.x = element_text(size = 25, face = "bold")  
  )



