# =====================================================================
# 6. LinDA Differential Abundance Analysis: Object Preprocessing
# =====================================================================

# ---------------------------------------------------------------------
# Step 1. Agglomeration, Renaming & Pre-processing
# ---------------------------------------------------------------------

if (!isTRUE(linda_params$skip_preprocessing)) {
  
  # 1a. Select base phyloseq object
  ps_linda <- linda_params$ps_input
  
  # 1b. Subsample input object if enabled (for pipeline testing)
  if (isTRUE(linda_params$subsample_input) && !is.null(linda_params$subsample_n)) {
    n_keep <- min(linda_params$subsample_n, nsamples(ps_linda))
    message("Subsampling input object to first ", n_keep, " samples for testing.")
    sub_names <- sample_names(ps_linda)[seq_len(n_keep)]
    ps_linda <- prune_samples(sub_names, ps_linda)
  }
  
  # 1c. Agglomerate to target taxonomic level
  target_level_linda <- linda_params$tax_level
  message("Agglomerating taxa to level: ", target_level_linda)
  ps_linda_tax <- tax_glom(ps_linda, taxrank = target_level_linda)
  
  # 1d. Make taxonomic names unique at the target level
  ps_linda_renamed <- renameTaxa(
    ps_linda_tax, 
    pat = "<name>", 
    substPat = "<subst_name>_<subst_R>",
    numDupli = target_level_linda,
    numDupliPat  = "<name>_<num>",  # Explicit separator prevents regex character-eating
    numUnclassPat = "<name>_<num>",
  )
  
  # 1e. Cache the processed object for direct loading in future runs
  if (!is.null(linda_params$cache_path)) {
    dir.create(dirname(linda_params$cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(ps_linda_renamed, file = linda_params$cache_path)
    message("Saved pre-processed object to: ", linda_params$cache_path)
  }
  
} else {
  # Load cached object directly if pre-processing step is skipped
  message("Skipping pre-processing. Loading cached object from: ", linda_params$cache_path)
  ps_linda_renamed <- readRDS(linda_params$cache_path)
}


# ---------------------------------------------------------------------
# Step 2. Dynamic Two-Layer Sample Subsetting
# ---------------------------------------------------------------------

ps_filtered <- ps_linda_renamed

# Layer 1 Filtering
if (!is.null(linda_params$filter1_var) && length(linda_params$filter1_val) > 0) {
  message("Applying Layer 1 filter: ", linda_params$filter1_var, " in [", paste(linda_params$filter1_val, collapse = ", "), "]")
  
  metadata <- sample_data(ps_filtered)
  keep_samples <- metadata[[linda_params$filter1_var]] %in% linda_params$filter1_val
  ps_filtered <- prune_samples(keep_samples, ps_filtered)
}

# Layer 2 Filtering (Optional)
if (!is.null(linda_params$filter2_var) && length(linda_params$filter2_val) > 0) {
  message("Applying Layer 2 filter: ", linda_params$filter2_var, " in [", paste(linda_params$filter2_val, collapse = ", "), "]")
  
  metadata <- sample_data(ps_filtered)
  keep_samples <- metadata[[linda_params$filter2_var]] %in% linda_params$filter2_val
  ps_filtered <- prune_samples(keep_samples, ps_filtered)
}

# Remove taxa that are no longer present in any of the subsetted samples
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)


# ---------------------------------------------------------------------
# Step 3. Prevalence Taxonomic Filtering
# ---------------------------------------------------------------------

if (!is.null(linda_params$min_prevalence) && linda_params$min_prevalence > 0) {
  
  n_samples_current <- nsamples(ps_filtered)
  prevalence_cutoff <- linda_params$min_prevalence * n_samples_current
  
  message(sprintf(
    "Prevalence threshold: present in >= %.1f%% of subset samples (%d out of %d samples)", 
    linda_params$min_prevalence * 100, 
    ceiling(prevalence_cutoff), 
    n_samples_current
  ))
  
  # Retain taxa meeting the minimum sample prevalence requirement
  core_taxa <- taxa_names(filter_taxa(ps_filtered, function(x) sum(x > 0) >= prevalence_cutoff, TRUE))
  ps_taxa_filtered <- prune_taxa(core_taxa, ps_filtered)
  
} else {
  ps_taxa_filtered <- ps_filtered
}

rm(ps_filtered)


# ---------------------------------------------------------------------
# Step 4. Metadata Cleaning & Formatting for LinDA
# ---------------------------------------------------------------------

metadata_correct_format <- data.frame(sample_data(ps_taxa_filtered))

categ_env_vars <- linda_params$categ_env_variables
numeric_env_vars <- linda_params$numeric_env_variables

# 4a. Automated Categorical / Factor Conversion
for (cat_var in categ_env_vars) {
  if (cat_var %in% colnames(metadata_correct_format)) {
    metadata_correct_format[[cat_var]] <- factor(as.character(metadata_correct_format[[cat_var]]))
    message(paste("Successfully converted to factor:", cat_var))
  } else {
    warning(paste("Configured categorical variable not found in dataset columns:", cat_var))
  }
}

# 4b. Automated Numerical Conversion Layer
for (num_var in numeric_env_vars) {
  if (num_var %in% colnames(metadata_correct_format)) {
    clean_chars <- trimws(as.character(metadata_correct_format[[num_var]]))
    
    suppressWarnings({
      numeric_vector <- as.numeric(clean_chars)
    })
    
    metadata_correct_format[[num_var]] <- numeric_vector
    message(paste("Successfully converted to numeric metric:", num_var))
    
    raw_na_count <- sum(is.na(clean_chars))
    new_na_count <- sum(is.na(numeric_vector))
    if (new_na_count > raw_na_count) {
      warning(paste0(
        "Variable '", num_var, "' contained non-numeric text strings (e.g., 'ND', '<LOD', or blanks) ",
        "which have been automatically forced to NA values."
      ))
    }
  } else {
    warning(paste("Configured numeric variable not found in dataset columns:", num_var))
  }
}

# 4c. Retain only targeted metadata variables
all_target_cols <- unique(c(categ_env_vars, numeric_env_vars))
valid_target_cols <- intersect(all_target_cols, colnames(metadata_correct_format))

linda_metadata <- metadata_correct_format[, valid_target_cols, drop = FALSE]

# ---------------------------------------------------------------------
# Step 5. Reconstruct Final Phyloseq Object for LinDA
# ---------------------------------------------------------------------

ps_for_linda_analysis <- ps_taxa_filtered
sample_data(ps_for_linda_analysis) <- sample_data(linda_metadata)

# Assign directly into global workspace
assign("ps_for_linda_analysis", ps_for_linda_analysis, envir = .GlobalEnv)

# Optional Cache Output
if (isTRUE(linda_params$save_final_rds) && !is.null(linda_params$final_rds_path)) {
  dir.create(dirname(linda_params$final_rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ps_for_linda_analysis, file = linda_params$final_rds_path)
  message("Saved final LinDA object to: ", linda_params$final_rds_path)
}

message("LinDA pre-processing complete. Output object ready: 'ps_for_linda_analysis'")