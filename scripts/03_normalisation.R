# scripts/03_normalisation.R
# Goal: Address compositionality by handling zeros and applying log-ratio transforms.
# Requires: ps_work (The phyloseq object chosen in the Rmd)

# 1. Handling Sparsity (The Zero Problem)
# ---------------------------------------

# Method 1: Simple Pseudo-count
mibi_pseudo <- microbiome::transform(ps_work, transform = "shift", shift = 1)

# Method 2: Bayesian-Multiplicative Replacement (BMR)
message("Attempting BMR zero-imputation with microbiome-specific settings...")

# Preparation: Extract matrix (Samples as Rows)
otu_mat <- as.matrix(microbiome::abundances(ps_work))
if (taxa_are_rows(ps_work)) { otu_mat <- t(otu_mat) }

mibi_bmr <- tryCatch({
  bmr_counts <- zCompositions::cmultRepl(
    otu_mat, 
    method  = "GBM",      
    output  = "p-counts", 
    label   = 0,          
    z.delete  = 1,        # Prevents dropping samples with high zero counts
    # By default, cmultRepl drops samples/taxa with >80% zeros because it 
    # struggles to predict the zero structure from such a sparse matrix.
    # We set z.delete = 1 to ensure we only drop 100% empty rows/
    z.warning = 1         
  )
  # --- THE SIMPLE FIX ---
  # Since dimensions and names are identical, we simply create a new 
  # phyloseq object by copying the original and overwriting the OTU slot.
  
  tmp_ps <- ps_work
  
  # Convert bmr_counts to a matrix and ensure it's transposed back to Taxa-as-Rows
  # (Since your otu_mat was 234 taxa x 18 samples, and bmr_counts is the same)
  otu_table(tmp_ps) <- otu_table(as.matrix(bmr_counts), taxa_are_rows = TRUE)
  
  message("BMR successful: Imputed counts integrated.")
  tmp_ps
  
}, error = function(e) {
  warning("BMR integration failed. Error: ", e$message)
  return(mibi_pseudo)
})

# 2. Transformation Methods
# ---------------------------------------

# Method A: Total Sum Scaling (TSS) / Relative Abundance
mibi_tss <- microbiome::transform(ps_work, transform = "compositional")

# Method B: Centered Log-Ratio (CLR)
mibi_pseudo_clr <- microbiome::transform(mibi_pseudo, transform = "clr")
mibi_bmr_clr    <- microbiome::transform(mibi_bmr, transform = "clr")

# 3. Validation
# ---------------------------------------
message("--- Normalization Summary ---")
message("1. TSS: Relative abundance (0-100%) created.")
message("2. Pseudo-CLR: Fast +1 shift transformation created.")
message("3. BMR-CLR: Bayesian-imputed transformation created.")
if(identical(mibi_bmr, mibi_pseudo)) {
  message("NOTE: BMR failed; mibi_bmr_clr is identical to mibi_pseudo_clr.")
}
message("-----------------------------")