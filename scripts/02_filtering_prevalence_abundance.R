# 02_filtering.R
# This script filters the phyloseq object based on prevalence and/or abundance

# 1. Define the filtering functions
# ---------------------------------

# Prevalence: ASV must be present in at least X samples
# Abundance: ASV must have at least Y reads in those samples

if (prev_filter == 1 && abund_filter == 0) {
  # ONLY PREVALENCE
  # Function: Sum of samples where count > 0 must be >= threshold
  mibi_phylo_object_filt <- filter_taxa(mibi_phylo_object, function(x) sum(x > 0) >= prev_filter_threshold, prune = TRUE)
  message("Filtering applied: Prevalence only.")
  
} else if (prev_filter == 0 && abund_filter == 1) {
  # ONLY ABUNDANCE
  # Function: ASV must have at least 'abund_filter_threshold' total reads across all samples
  mibi_phylo_object_filt <- filter_taxa(mibi_phylo_object, function(x) sum(x) >= abund_filter_threshold, prune = TRUE)
  message("Filtering applied: Abundance only.")
  
} else if (prev_filter == 1 && abund_filter == 1) {
  # BOTH (Strict): Must meet abundance AND prevalence
  # Function: Number of samples where count >= abundance_threshold must be >= prev_threshold
  mibi_phylo_object_filt <- filter_taxa(mibi_phylo_object, function(x) sum(x >= abund_filter_threshold) >= prev_filter_threshold, prune = TRUE)
  message("Filtering applied: Combined Prevalence and Abundance.")
  
} else {
  # NO FILTERING
  mibi_phylo_object_filt <- mibi_phylo_object
  message("No filtering applied (Check switches).")
}

# 2. Print result for the report
ntaxa_before <- ntaxa(mibi_phylo_object)
ntaxa_after <- ntaxa(mibi_phylo_object_filt)
cat("Filtered", ntaxa_before - ntaxa_after, "ASVs. Remaining:", ntaxa_after, "\n")