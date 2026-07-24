# 02_filtering_variance.R
# Optimized for Compositional Microbiome Data using IQR

# Ensure a starting object exists
if (!exists("mibi_phylo_object_filt")) {
  if (exists("mibi_phylo_object")) {
    mibi_phylo_object_filt <- mibi_phylo_object
  } else {
    stop("Error: No phyloseq object found to filter.")
  }
}

if (variance_filter_switch == 1) {
  
  # 1. Extract matrix (taxa are rows)
  otu <- as.matrix(microbiome::abundances(mibi_phylo_object_filt))
  
  # 2. Calculate IQR (The robust choice for compositional sparsity)
  # 
  taxon_variance <- apply(otu, 1, IQR)
  
  # 3. Identify Threshold
  thresh_val <- quantile(taxon_variance, probs = variance_cutoff, na.rm = TRUE)
  
  # 4. Filter Logic
  # Safeguard: If thresh_val is 0, we strictly remove everything with 0 variance.
  # Otherwise, we use the standard percentile.
  if (thresh_val == 0) {
    taxa_to_keep <- names(taxon_variance[taxon_variance > 0])
    note_msg <- " (Threshold was 0; dropped all absolute zero-variance taxa)"
  } else {
    taxa_to_keep <- names(taxon_variance[taxon_variance >= thresh_val])
    note_msg <- ""
  }
  
  # 5. Execute Pruning
  ntaxa_before_var <- ntaxa(mibi_phylo_object_filt)
  mibi_phylo_object_filt <- prune_taxa(taxa_to_keep, mibi_phylo_object_filt)
  ntaxa_after_var <- ntaxa(mibi_phylo_object_filt)
  
  # 6. Reporting
  message(paste0("Filtering applied: Low variance (Interquartile Range).", note_msg))
  cat("Variance Filter: Excluded the bottom ", variance_cutoff * 100, 
      "% taxa. Removed: ", ntaxa_before_var - ntaxa_after_var, " ASVs.\n", sep="")
  
} else {
  message("No variance filtering applied (Switch is Off).")
}

cat("Final taxa count:", ntaxa(mibi_phylo_object_filt), "\n")