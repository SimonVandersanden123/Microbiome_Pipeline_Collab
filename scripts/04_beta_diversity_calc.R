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

# --- 4. Statistical Testing (PERMANOVA with Streamlined Factor Processing) ---
# ----------------------------------------------------------------------------
dist_matrix <- phyloseq::distance(ps_beta_input, method = dist_method)
metadata    <- as(sample_data(ps_beta_input), "data.frame")

# Combine all targeted variables for a complete case check
all_target_vars <- c(numeric_env_variables, categ_env_variables)

# 1. Create a temporary matrix strictly to find complete cases (no NAs)
# We force numeric columns here so that text values like "missing" turn to NA and get dropped
# Here we need to watch out for . instead of , in the metadata
temp_check <- metadata %>%
  dplyr::select(dplyr::all_of(all_target_vars)) %>%
  mutate(across(dplyr::all_of(numeric_env_variables), function(x) as.numeric(as.character(x))))

complete_indices    <- which(complete.cases(temp_check))
complete_sample_ids <- rownames(metadata)[complete_indices]

# 2. Subset BOTH your final metadata and distance matrix to match perfectly
metadata_complete   <- metadata[complete_sample_ids, , drop = FALSE]
dist_matrix_complete <- as.dist(as.matrix(dist_matrix)[complete_sample_ids, complete_sample_ids])

# 3. Apply clean data types to your final metadata object (THE ONLY LOOP YOU NEED)
metadata_complete[numeric_env_variables] <- lapply(metadata_complete[numeric_env_variables], function(x) as.numeric(as.character(x)))

# =====================================================================
# AUTOMATED FACTOR CONVERSION LAYER (Generalized)
# =====================================================================
for (cat_var in categ_env_variables) {
  if (cat_var %in% colnames(metadata_complete)) {
    
    # Clean up column to a character format first, then cast to a standard factor
    # This ensures R handles any messy mixed-type inputs or text labels cleanly
    metadata_complete[[cat_var]] <- factor(as.character(metadata_complete[[cat_var]]))
    
    message(paste("Successfully converted to factor:", cat_var))
    
  } else {
    warning(paste("Configured categorical variable not found in dataset columns:", cat_var))
  }
}
# =====================================================================
# =====================================================================
# 4. Formulate formula dynamically and run PERMANOVA

formula_beta       <- reformulate(all_target_vars, response = "dist_matrix_complete")
message("Running marginal PERMANOVA using streamlined numeric and factor matrices...")
permanova_marginal <- adonis2(formula_beta, data = metadata_complete, by = "margin", permutations = 999)
#this is just to test all varriables seperatly, while taking into account the other variables.
# print(permanova_marginal) # just give the permanova results, in the order in which the variables were put in
# Print the Permanova results ordered from highest R2 value
# 1. Convert the PERMANOVA object to a clean, sortable dataframe
permanova_sorted <- as.data.frame(permanova_marginal) %>%
  tibble::rownames_to_column("Variable") %>%
  # 2. Separate your actual variables from the Residual and Total rows
  filter(!Variable %in% c("Residual", "Total")) %>%
  # 3. Sort by R2 in descending order (highest on top)
  arrange(desc(R2))

# 4. Bind the Residual and Total rows back to the bottom so the math remains intact
structural_rows <- as.data.frame(permanova_marginal) %>%
  tibble::rownames_to_column("Variable") %>%
  filter(Variable %in% c("Residual", "Total"))

permanova_final_table <- bind_rows(permanova_sorted, structural_rows)
beta_p_val <- permanova_final_table$`Pr(>F)`[1]
# 5. Print your beautifully ordered table!
print(permanova_final_table, row.names = FALSE)

# 2. Conditionally run the sequential PERMANOVA if the config variable exists
# We check if the variable exists and is not empty
if (!is.null(config$Beta_Diversity$specified_permanova_formula) && 
    length(config$Beta_Diversity$specified_permanova_formula) > 0) {
  
  message("Specific formula detected. Running sequential PERMANOVA...")

  # Retrieve the formula terms from config
  formula_terms <- config$Beta_Diversity$specified_permanova_formula
  # Create the formula object
  # Note: reformulate expects a character vector of terms
  formula_beta_terms <- reformulate(formula_terms, response = "dist_matrix_complete")

  # Run adonis2 with by = "terms"
  permanova_terms <- adonis2(formula_beta_terms, data = metadata_complete, by = "terms", permutations = 999)
  
  #Now clean up the results before printing the output
  # 1. Convert the PERMANOVA object to a clean, sortable dataframe
  permanova_sorted_terms <- as.data.frame(permanova_terms) %>%
    tibble::rownames_to_column("Variable") %>%
    # 2. Separate your actual variables from the Residual and Total rows
    filter(!Variable %in% c("Residual", "Total")) %>%
    # 3. Sort by R2 in descending order (highest on top)
    arrange(desc(R2))
  
  # 4. Bind the Residual and Total rows back to the bottom so the math remains intact
  structural_rows_terms <- as.data.frame(permanova_terms) %>%
    tibble::rownames_to_column("Variable") %>%
    filter(Variable %in% c("Residual", "Total"))
  # 5. Print your beautifully ordered table!
  print(permanova_terms_final_table, row.names = FALSE)
  
  # Extract p-value for the plot title later
  beta_p_val <- permanova_sorted_terms$`Pr(>F)`[1]
  
  permanova_terms_final_table <- bind_rows(permanova_sorted_terms, structural_rows_terms)
} else {
  message("No specific PERMANOVA formula provided in config. Skipping sequential analysis.")
}


