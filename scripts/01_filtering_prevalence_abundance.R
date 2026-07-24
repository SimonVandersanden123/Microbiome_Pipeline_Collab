# scripts/02_filtering_prevalence_abundance.R
# Requires: mibi_phylo_object, prev_filter_threshold, abund_filter_threshold...

# 1. Validation check: ensure variables exist
if (!exists("prev_filter_threshold")) stop("Please define prev_filter_threshold in the Rmd.")

# 2. Execution logic
if (prev_filter == 1 && abund_filter == 0) {
  mibi_phylo_object_filt <- filter_taxa(mibi_phylo_object, function(x) sum(x > 0) >= prev_filter_threshold, prune = TRUE)
  filter_msg <- paste("Prevalence only (min samples:", prev_filter_threshold, ")")
  
} else if (prev_filter == 0 && abund_filter == 1) {
  mibi_phylo_object_filt <- filter_taxa(mibi_phylo_object, function(x) sum(x) >= abund_filter_threshold, prune = TRUE)
  filter_msg <- paste("Abundance only (min total reads:", abund_filter_threshold, ")")
  
} else if (prev_filter == 1 && abund_filter == 1) {
  # Strict Combined Filter
  mibi_phylo_object_filt <- filter_taxa(mibi_phylo_object, function(x) sum(x >= abund_filter_threshold) >= prev_filter_threshold, prune = TRUE)
  filter_msg <- paste("Combined (min", abund_filter_threshold, "reads in", prev_filter_threshold, "samples)")
  
} else {
  mibi_phylo_object_filt <- mibi_phylo_object
  filter_msg <- "None applied."
}

# 3. Summary Stats for the User
ntaxa_before <- ntaxa(mibi_phylo_object)
ntaxa_after <- ntaxa(mibi_phylo_object_filt)
percent_kept <- round((ntaxa_after / ntaxa_before) * 100, 2)

message("--- Filtering Summary ---")
message("Criteria: ", filter_msg)
message("Taxa before: ", ntaxa_before)
message("Taxa after:  ", ntaxa_after)
message("Percent kept: ", percent_kept, "%")
message("-------------------------")
