
#Handling the "Zero Problem" and Compositional Transformations

library(microbiome) # Excellent for CLR and TSS functions

# 1. Handling Sparsity (The Zero Problem)
# ---------------------------------------
# For CLR, we must replace zeros. Here we add a pseudo-count of 1.
# This shifts the data into a space where log(x) is defined.

mibi_pseudo <- transform(mibi_phylo_object_filt, transform = "shift", shift = 1)

#Add a part when proper pseudocounting, using Bayesian-Multiplicative Replacement

# 2. Transformation Methods
# ---------------------------------------

# Method A: Total Sum Scaling (TSS) / Relative Abundance
# Logic: count(i) / sum(all counts in sample)
mibi_tss <- transform(mibi_phylo_object_filt, transform = "compositional")

# Method B: Centered Log-Ratio (CLR)
# Logic: log(x_i / geometric_mean(x))
# This moves data from the "simplex" (bounded 0-1) to Euclidean space.
mibi_clr <- transform(mibi_pseudo, transform = "clr")

# 3. Validation
# ---------------------------------------
message("Normalization complete:")
message("- TSS object created (for visualization)")
message("- CLR object created (for linear modeling/correlation)")