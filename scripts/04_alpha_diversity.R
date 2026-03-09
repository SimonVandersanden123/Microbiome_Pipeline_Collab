# This script calculates different alpha diversity metrics
# 03_alpha_diversity.R
library(microbiome)
library(picante)

# 1. Master Calculation
# ---------------------
# Calculate every available metric in one go to ensure consistency
tab_all <- microbiome::alpha(alfa_diversity_phyloseq_object, index = "all")

# 2. Category Subsetting Logic
# ----------------------------
# We map your user-defined vectors to the column names in the master table
# Note: microbiome::alpha uses 'diversity_shannon' for shannon, etc.

# Define the mapping (adjusting for microbiome::alpha naming conventions)
rich_cols <- paste0("observed") # or "chao1"
div_cols  <- paste0("diversity_", selected_information)
even_cols <- paste0("evenness_", selected_evenness)
dom_cols  <- paste0("dominance_", tolower(selected_dominance))

# Combine all requested columns
all_requested_cols <- c(selected_richness, div_cols, even_cols, dom_cols)

# Subset the master table
alpha_final <- tab_all[, colnames(tab_all) %in% all_requested_cols, drop = FALSE]

# 4. Final Merge with Metadata
# ----------------------------
alpha_output <- cbind(sample_data(alfa_diversity_phyloseq_object), alpha_final)

message("Alpha diversity calculation complete. Table generated for selected traits.")