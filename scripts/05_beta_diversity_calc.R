# scripts/05_beta_diversity_calc.R

# --- 1. Object Selection Logic ---
if (beta_metric == "aitchison") {
  # Aitchison = Euclidean distance on CLR data
  ps_beta_input <- if(clr_variant == "BMR") mibi_bmr_clr else mibi_pseudo_clr
  dist_method <- "euclidean"
  ord_method <- "RDA" # RDA with no constraints = PCA
  msg_text <- paste0("Aitchison Distance (PCA) using ", clr_variant, " imputation")
} else {
  # Standard distances on Relative Abundance (TSS)
  ps_beta_input <- mibi_tss
  dist_method <- beta_metric
  ord_method <- "PCoA"
  msg_text <- paste0(beta_metric, " Distance (PCoA) on TSS data")
}

# --- 2. Safety Check (Phylogeny) ---
if (grepl("unifrac", beta_metric, ignore.case = TRUE) & is.null(phy_tree(ps_beta_input))) {
  stop("Error: UniFrac requires a phylogenetic tree. None found in object.")
}
# --- 2. Safety Check (Phylogeny) ---
if (grepl("wunifrac", beta_metric, ignore.case = TRUE) & is.null(phy_tree(ps_beta_input))) {
  stop("Error: wuniFrac requires a phylogenetic tree. None found in object.")
}

# --- 3. Run Ordination ---
message(paste("Calculating:", msg_text))
ord_beta <- ordinate(ps_beta_input, method = ord_method, distance = dist_method)

# --- 4. Statistical Testing (PERMANOVA) ---
dist_matrix <- phyloseq::distance(ps_beta_input, method = dist_method)
metadata <- as(sample_data(ps_beta_input), "data.frame")

formula_beta <- reformulate(color_var, response = "dist_matrix")
permanova_res <- adonis2(formula_beta, data = metadata)

# Extract p-value for the plot title later
beta_p_val <- permanova_res$`Pr(>F)`[1]