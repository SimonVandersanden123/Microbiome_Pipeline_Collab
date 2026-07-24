# 08_barchart_specific_group.R

#' Get best available taxonomic name (Fallback logic)
get_clean_name <- function(tax_row) {
  levels <- c("genus", "family", "order", "class", "phylum")
  # Prefixes help identify what level the fallback stopped at
  prefixes <- c("", "f_", "o_", "c_", "p_")
  
  for (i in seq_along(levels)) {
    val <- tax_row[[levels[i]]]
    if (!is.na(val) && !str_detect(tolower(val), "unclassified|incertae") && val != "") {
      return(paste0(prefixes[i], val))
    }
  }
  return("Unclassified")
}

#' Subset and Prepare Specific Taxa with Shaded Hierarchy
#' @param ps Phyloseq object
#' @param rank_name The rank to filter by (e.g., "phylum")
#' @param rank_value The specific taxon to keep (e.g., "Bacteroidota")
#' @param high_level The grouping level for colors (e.g., "family")
#' @param low_level The detailed level for labels (e.g., "genus")
#' @param facet_var Metadata column for grouping
#' @param threshold % Abundance threshold within this subset
prepare_specific_abundance_data <- function(ps, rank_name, rank_value, 
                                            high_level = "family", low_level = "genus", 
                                            facet_var = "Sample", threshold = 1) {
  
  # 1. Subset to the specific group and transform to %
  tax <- as.data.frame(tax_table(ps))
  keep_otus <- rownames(tax)[tax[[rank_name]] == rank_value]
  # Handle cases where the taxon name doesn't exist to avoid a crash
  if(length(keep_otus) == 0) stop(paste("Taxon", rank_value, "not found in rank", rank_name))
  
  ps_sub <- prune_taxa(keep_otus, ps)
  # 2. Transform to % (Relative to this subset)
  ps_rel <- transform_sample_counts(ps_sub, function(x) x / sum(x) * 100)
  df <- psmelt(ps_rel)
  
  # 2. Average by group if specified
  if (facet_var != "Sample") {
    df <- df %>%
      dplyr::group_by(OTU, !!sym(facet_var), !!sym(high_level), !!sym(low_level)) %>%
      dplyr::summarise(Abundance = mean(Abundance), .groups = "drop")
  }
  
  # 3. Create Hierarchical Labels and handle "Others"
  df <- df %>%
    dplyr::mutate(
      Level_High = as.character(!!sym(high_level)),
      Level_Low  = as.character(!!sym(low_level)),
      # If name is missing, use the "Clean Name" fallback logic
      Plot_High = ifelse(Abundance >= threshold, Level_High, "Others"),
      Plot_Low  = ifelse(Abundance >= threshold, Level_Low, "Minor Taxa"),
      Hierarchical_Label = ifelse(Plot_High == "Others", "Others", paste(Plot_High, Plot_Low, sep = " - "))
    )
  
  # 4. Collapse OTUs into solid blocks
  df <- df %>%
    dplyr::group_by(!!sym(facet_var), Plot_High, Hierarchical_Label) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # 5. Order "Others" to the bottom
  all_labels <- unique(df$Hierarchical_Label)
  major_taxa <- sort(setdiff(all_labels, "Others"))
  target_levels <- c(major_taxa, "Others")
  
  df$Hierarchical_Label <- factor(df$Hierarchical_Label, levels = target_levels)
  
  return(df)
}