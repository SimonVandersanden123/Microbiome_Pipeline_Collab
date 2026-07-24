# This script calculates different alpha diversity metrics
# 04_alpha_diversity.R
# 1. Master Calculation
tab_all <- microbiome::alpha(alfa_diversity_phyloseq_object, index = "all")

# 2. Category Mapping Logic
# We need to match the specific prefixes used by the microbiome package
requested_metrics <- c(
  selected_richness, 
  paste0("diversity_", selected_information),
  paste0("evenness_", selected_evenness),
  paste0("dominance_", tolower(selected_dominance))
)

# 3. Subset and clean names
# Only keep columns that actually exist in the output
alpha_final <- tab_all[, colnames(tab_all) %in% requested_metrics, drop = FALSE]

# Simplify names for the plot (removing 'diversity_' prefix etc. for cleaner labels)
colnames(alpha_final) <- gsub("diversity_|evenness_|dominance_", "", colnames(alpha_final))
plot_metrics_all <- colnames(alpha_final) # Save for the plotting script

# 4. Merge with Metadata
alpha_output <- cbind(as(sample_data(alfa_diversity_phyloseq_object), "data.frame"), alpha_final)

message("Alpha diversity calculation complete.")