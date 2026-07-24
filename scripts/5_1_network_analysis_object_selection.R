# =====================================================================
# 5. Microbial Network Analysis: Object Preparation & Filtering
# =====================================================================

library(NetCoMi)
library(phyloseq)

# ---------------------------------------------------------------------
# Step 1. Pre-processing & Taxonomic Agglomeration
# ---------------------------------------------------------------------

if (!isTRUE(net_params$skip_preprocessing)) {
  
  # 1a. Select base phyloseq object
  ps_net <- net_params$ps_input
  
  # 1b. Subsample input object if enabled (e.g. for pipeline testing)
  if (isTRUE(net_params$subsample_input) && !is.null(net_params$subsample_n)) {
    n_keep <- min(net_params$subsample_n, nsamples(ps_net))
    message("Subsampling input object to first ", n_keep, " samples.")
    sub_names <- sample_names(ps_net)[seq_len(n_keep)]
    ps_net <- prune_samples(sub_names, ps_net)
  }
  
  # 1c. Agglomerate to target taxonomic level
  message("Agglomerating taxa to level: ", net_params$tax_level)
  ps_net_tax <- tax_glom(ps_net, taxrank = net_params$tax_level)
  
  # 1d. Make taxonomic names unique at the target level for NetCoMi
  ps_net_renamed <- renameTaxa(
    ps_net_tax, 
    pat = "<name>", 
    substPat = "<name>_<subst_name>(<subst_R>)",
    numDupli = net_params$tax_level,
    numDupliPat  = "<name>_<num>",  # Explicit separator prevents regex character-eating
    numUnclassPat = "<name>_<num>",
  )
  
  # 1e. Cache the processed object for direct loading in future runs
  if (!is.null(net_params$cache_path)) {
    dir.create(dirname(net_params$cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(ps_net_renamed, file = net_params$cache_path)
    message("Saved pre-processed object to: ", net_params$cache_path)
  }
  
} else {
  # Load cached object directly if pre-processing step is skipped
  message("Skipping pre-processing. Loading cached object from: ", net_params$cache_path)
  ps_net_renamed <- readRDS(net_params$cache_path)
}


# ---------------------------------------------------------------------
# Step 2. Dynamic Two-Layer Sample Subsetting
# ---------------------------------------------------------------------

ps_filtered <- ps_net_renamed

# Layer 1 Filtering
if (!is.null(net_params$filter1_var) && length(net_params$filter1_val) > 0) {
  message("Applying Layer 1 filter: ", net_params$filter1_var, " in [", paste(net_params$filter1_val, collapse = ", "), "]")
  
  metadata <- sample_data(ps_filtered)
  keep_samples <- metadata[[net_params$filter1_var]] %in% net_params$filter1_val
  ps_filtered <- prune_samples(keep_samples, ps_filtered)
}

# Layer 2 Filtering (Optional)
if (!is.null(net_params$filter2_var) && length(net_params$filter2_val) > 0) {
  message("Applying Layer 2 filter: ", net_params$filter2_var, " in [", paste(net_params$filter2_val, collapse = ", "), "]")
  
  metadata <- sample_data(ps_filtered)
  keep_samples <- metadata[[net_params$filter2_var]] %in% net_params$filter2_val
  ps_filtered <- prune_samples(keep_samples, ps_filtered)
}

# Remove taxa that are no longer present in any of the subsetted samples
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)


# ---------------------------------------------------------------------
# Step 3. Prevalence Taxonomic Filtering
# ---------------------------------------------------------------------

if (!is.null(net_params$min_prevalence) && net_params$min_prevalence > 0) {
  
  n_samples_current <- nsamples(ps_filtered)
  prevalence_cutoff <- net_params$min_prevalence * n_samples_current
  
  message(sprintf(
    "Prevalence threshold: present in >= %.1f%% of subset samples (%d out of %d samples)", 
    net_params$min_prevalence * 100, 
    ceiling(prevalence_cutoff), 
    n_samples_current
  ))
  
  # Retain taxa meeting the minimum sample prevalence requirement
  core_taxa <- taxa_names(filter_taxa(ps_filtered, function(x) sum(x > 0) >= prevalence_cutoff, TRUE))
  ps_object_for_network_analysis <- prune_taxa(core_taxa, ps_filtered)
  
} else {
  ps_object_for_network_analysis <- ps_filtered
}

# Clean up execution memory
rm(ps_filtered)

message("Network analysis object successfully created: 'ps_object_for_network_analysis'")


# =====================================================================
# Original Raw Reference Script (v1.0)
# =====================================================================
# library(NetCoMi)
# # =====================================================================
# # 1. Selection and filtering of the phyloseq object to use for network selection
# # =====================================================================
# # First: Select the phyloseq object choose the original filtered for interest data table with no prior normalisation.
# ps_network_analysis <- ps_work
# first_50_samples <- sample_names(ps_network_analysis)[1:162]
# ps_debug <- subset_samples(ps_network_analysis, sample_names(ps_network_analysis) %in% first_50_samples)
# 
# otu_debug <- as(otu_table(ps_debug), "matrix")
# ps_debug_genus <- tax_glom(ps_debug, taxrank = "genus")
# 
# ps_debug_genus_renamed <- renameTaxa(ps_debug_genus, 
#                                      pat = "<name>", 
#                                      substPat = "<name>_<subst_name>(<subst_R>)",
#                                      numDupli = "genus")
# saveRDS(ps_debug_genus_renamed, file = "Objects_Generated_In_Pipeline/Network_Analysis/ps_network_analysis_Genus_renamed_wetland_samples.rds")
# 
# ps_debug_genus_renamed <- readRDS("Objects_Generated_In_Pipeline/Network_Analysis/ps_network_analysis_Genus_renamed_wetland_samples.rds")
# 
# # ====================================================================
# # 2. Subsetting and filtering
# # ====================================================================
# ps_debug_genus_renamed_selection <- subset_samples(ps_debug_genus_renamed, timepoint %in% c("T0","T1","T2","T3"))
# ps_debug_genus_renamed_selection_2layer_filtering <- subset_samples(ps_debug_genus_renamed_selection, CW_number %in% c("CW_2"))
# ps_debug_genus_renamed_selection <- ps_debug_genus_renamed_selection_2layer_filtering
# 
# # =====================================================================
# # Further filtering for the network analysis, prevalence 
# # =====================================================================
# prevalence_threshold <- 0.25 * nsamples(ps_debug_genus_renamed_selection)
# core_genera <- taxa_names(filter_taxa(ps_debug_genus_renamed_selection, function(x) sum(x > 0) >= prevalence_threshold, TRUE))
# ps_prevalence_filtered <- prune_taxa(core_genera, ps_debug_genus_renamed_selection)
# ps_object_for_network_analysis <- ps_prevalence_filtered