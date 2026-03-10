# 05_beta_diversity_calc.R
# Universal Beta Diversity Calculation (Aitchison, Bray, Jaccard, UniFrac)

library(phyloseq)
library(microbiome)
library(vegan)

# --- 1. Object Selection Logic ---
if (beta_metric == "aitchison") {
  # Use the CLR-transformed objects for Aitchison (Euclidean on CLR)
  ps_work <- if(clr_variant == "BMR") mibi_bmr_clr else mibi_pseudo_clr
  dist_method <- "euclidean"
  ord_method <- "RDA" # RDA with no constraints = PCA
  msg_text <- paste0("Aitchison Distance (PCA) using ", clr_variant, " imputation")
} else {
  # Use Relative Abundance (TSS) for Jaccard, Bray, and UniFrac
  ps_work <- mibi_tss
  dist_method <- beta_metric
  ord_method <- "PCoA"
  msg_text <- paste0(beta_metric, " Distance (PCoA) on TSS data")
}

# --- 2. Check for Phylogenetic Tree (Safety for UniFrac) ---
if (grepl("unifrac", beta_metric, ignore.case = TRUE) & is.null(phy_tree(ps_work))) {
  stop("Error: You selected UniFrac, but the phyloseq object has no phylogenetic tree.")
}

# --- 3. Run Ordination ---
message(paste("Calculating:", msg_text))
ord_beta <- ordinate(ps_work, method = ord_method, distance = dist_method)

# --- 4. Statistical Testing (PERMANOVA) ---
# Calculate the distance matrix for adonis2
dist_matrix <- phyloseq::distance(ps_work, method = dist_method)
metadata <- as(sample_data(ps_work), "data.frame")

# Use reformulate to handle the user-defined color_var string
formula_beta <- reformulate(color_var, response = "dist_matrix")
permanova_res <- adonis2(formula_beta, data = metadata)

message("Calculation Complete. Use the plotting chunk to visualize.")