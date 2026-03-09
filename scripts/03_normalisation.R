
#Handling the "Zero Problem" and Compositional Transformations

library(microbiome) # Excellent for CLR and TSS functions
library(zCompositions)

# 1. Handling Sparsity (The Zero Problem)
# ---------------------------------------
# Method 1: Simple Pseudo-count
# Logic: Adds a constant (1) to all cells. Fast but can bias ratios.
mibi_pseudo <- transform(mibi_phylo_object_filt, transform = "shift", shift = 1)
######################################################################################
# Method 2: Bayesian-Multiplicative Replacement (BMR)
# Logic: Imputes zeros based on a Dirichlet prior and shifts non-zeros multiplicatively.
# Note: cmultRepl works on a samples-as-rows, taxa-as-columns matrix.
otu_table_matrix <- as.matrix(abundances(mibi_phylo_object_filt))
# Transpose because cmultRepl expects samples in rows
otu_mat_t <- t(otu_table_matrix)
# Perform BMR (Geometric Bayesian Multiplicative method)
# be carefull for singletons, it wont work with singletons
# output = "p-counts" returns pseudo-counts that preserve ratios
bmr_counts <- cmultRepl(otu_mat_t, method = "GBM", output = "p-counts")

# Re-integrate into a phyloseq object
mibi_bmr <- mibi_phylo_object_filt
otu_table(mibi_bmr) <- otu_table(t(bmr_counts), taxa_are_rows = TRUE)
######################################################################################
# Method 3: Matrix completion (MC)


# 2. Transformation Methods
# ---------------------------------------

# Method A: Total Sum Scaling (TSS) / Relative Abundance
# Logic: count(i) / sum(all counts in sample)
mibi_tss <- transform(mibi_phylo_object_filt, transform = "compositional")

# Method B: Centered Log-Ratio (CLR)
# Logic: log(x_i / geometric_mean(x))
# This moves data from the "simplex" (bounded 0-1) to Euclidean space.
mibi_pseudo_clr <- transform(mibi_pseudo, transform = "clr")

mibi_bmr_clr <- transform(mibi_bmr, transform = "clr")

mibi_mc_clr <- transform(mibi_mc, transform = "clr")
# 3. Validation
# ---------------------------------------
message("Normalization complete:")
message("- TSS object created (for visualization)")
message("- CLR object created (for linear modeling/correlation)")