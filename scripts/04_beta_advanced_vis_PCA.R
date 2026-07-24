# scripts/05_beta_advanced_vis_PCA.R

# =====================================================================
# 1. Coordinate & Variance Extraction
# =====================================================================
samp_coords <- as.data.frame(vegan::scores(ord_beta, display = "sites", choices = c(1, 2)))
taxa_coords <- as.data.frame(vegan::scores(ord_beta, display = "species", choices = c(1, 2)))

colnames(samp_coords) <- c("Dim1", "Dim2")
colnames(taxa_coords) <- c("Dim1", "Dim2")

# Calculate Variance Explained for Axis Labels
pca_summary <- summary(eigenvals(ord_beta))
pc1_var     <- round(pca_summary[2, 1] * 100, 1) # Proportion Explained for PC1
pc2_var     <- round(pca_summary[2, 2] * 100, 1) # Proportion Explained for PC2

# Merge sample coords with metadata for plotting
metadata_df          <- metadata_complete
samp_coords$SampleID <- rownames(samp_coords)
metadata_df$SampleID <- rownames(metadata_df)
ordination_df        <- merge(samp_coords, metadata_df, by = "SampleID")

# Determine max limit of sample points for auto-scaling
max_point_limit <- max(abs(c(ordination_df$Dim1, ordination_df$Dim2)))


# =====================================================================
# 2. Environmental Fitting & Centroids
# =====================================================================
# Select and order environmental variables to strictly match ordination samples
ord_sample_ids <- rownames(samp_coords) 
env_data       <- metadata_df[ord_sample_ids, c(numeric_env_variables, categ_env_variables), drop = FALSE]

# Fit environmental vectors and factors
enfit_env  <- vegan::envfit(ord_beta, env_data, permutations = 999, na.rm = TRUE)
env_coords <- as.data.frame(vegan::scores(enfit_env, "vectors")) 

if (nrow(env_coords) > 0) {
  env_coords$Variable <- rownames(env_coords)
  colnames(env_coords)[1:2] <- c("Dim1", "Dim2")
}

# Apply Conditional Filtering based on PERMANOVA P-values
if (isTRUE(filtering_env_variables_perm)) {
  permanova_df <- as.data.frame(permanova_marginal) %>% 
    tibble::rownames_to_column("Variable")
  
  sig_permanova_vars <- permanova_df %>%
    dplyr::filter(`Pr(>F)` < arrow_env_permanova_significance) %>%
    dplyr::pull(Variable)
  
  dropped_vars   <- setdiff(permanova_df$Variable, sig_permanova_vars)
  sig_env_arrows <- env_coords[env_coords$Variable %in% sig_permanova_vars, ]
  
  cat(sprintf(
    "\n[Environmental Vector Filter: ENABLED - p < %.3f]:\n - Total input variables: %d\n - Retained vectors      : %d (%s)\n - Dropped vectors       : %d (%s)\n\n",
    arrow_env_permanova_significance,
    nrow(permanova_df),
    length(sig_permanova_vars),
    if(length(sig_permanova_vars) > 0) paste(sig_permanova_vars, collapse = ", ") else "None",
    length(dropped_vars),
    if(length(dropped_vars) > 0) paste(dropped_vars, collapse = ", ") else "None"
  ))
} else {
  # Retain all variables if filtering is disabled
  sig_permanova_vars <- colnames(env_data)
  sig_env_arrows     <- env_coords
  
  cat(sprintf(
    "\n[Environmental Vector Filter: DISABLED]:\n - Retaining all %d environmental variables.\n\n",
    nrow(sig_env_arrows)
  ))
}  

# Auto-scale environmental arrows relative to the sample spread (~75% max extent)
if (nrow(sig_env_arrows) > 0) {
  scale_factor_env    <- (max_point_limit / max(abs(c(sig_env_arrows$Dim1, sig_env_arrows$Dim2)))) * 0.75
  sig_env_arrows$Dim1 <- sig_env_arrows$Dim1 * scale_factor_env
  sig_env_arrows$Dim2 <- sig_env_arrows$Dim2 * scale_factor_env
}

# Extract and Filter Categorical Centroids
env_centroids <- data.frame()

if (!is.null(enfit_env$factors)) {
  raw_centroids <- as.data.frame(vegan::scores(enfit_env, "factors"))
  
  if (nrow(raw_centroids) > 0) {
    raw_centroids$Level <- rownames(raw_centroids)
    colnames(raw_centroids)[1:2] <- c("Dim1", "Dim2")
    
    # Map factor levels back to parent variables
    raw_centroids$Variable <- NA
    for (var in sig_permanova_vars) {
      matches <- grepl(paste0("^", var), raw_centroids$Level)
      if (any(matches)) {
        raw_centroids$Variable[matches] <- var
      }
    }
    
    # Filter to globally significant variables
    sig_env_centroids <- raw_centroids %>% dplyr::filter(!is.na(Variable))
    
    if (nrow(sig_env_centroids) > 0) {
      # Strip variable prefix for display labels (e.g., "timepointPhase_1" -> "Phase_1")
      sig_env_centroids$Clean_Label <- sapply(1:nrow(sig_env_centroids), function(i) {
        gsub(paste0("^", sig_env_centroids$Variable[i]), "", sig_env_centroids$Level[i])
      })
      env_centroids <- sig_env_centroids
    }
  }
}


# =====================================================================
# 3. Taxa Loadings Logic & Filtering
# =====================================================================
tax_table_df             <- as.data.frame(tax_table(ps_beta_input))
taxa_coords$Taxon_Label  <- tax_table_df[[target_level]]
taxa_coords$Taxon_Label[is.na(taxa_coords$Taxon_Label)] <- rownames(taxa_coords)[is.na(taxa_coords$Taxon_Label)]

total_initial_taxa <- nrow(taxa_coords)

if (isTRUE(filter_taxon_loadings)) {
  message("Applying abundance and prevalence filtering to taxa arrows...")
  
  ps_for_filter  <- mibi_tss # Normalized TSS object
  otu_filter_mat <- as(phyloseq::otu_table(ps_for_filter), "matrix")
  if (!phyloseq::taxa_are_rows(ps_for_filter)) { otu_filter_mat <- t(otu_filter_mat) }
  
  asv_means       <- rowMeans(otu_filter_mat)
  asv_prevalences <- rowSums(otu_filter_mat > 0) / ncol(otu_filter_mat)
  
  n_samples                    <- nrow(metadata_complete)
  min_sample_count             <- config$Beta_Diversity$advanced$filter_taxon_prevalence_filter
  calculated_prevalence_filter <- min_sample_count / n_samples
  
  keep_asvs <- names(asv_means[asv_means >= filter_taxon_abundance_filter & 
                                 asv_prevalences >= calculated_prevalence_filter])
  
  taxa_coords_filtered <- taxa_coords[rownames(taxa_coords) %in% keep_asvs, ]
  
  taxa_kept    <- nrow(taxa_coords_filtered)
  taxa_removed <- total_initial_taxa - taxa_kept
  
  cat(sprintf("\n[Taxa Filter Report]:\n - Total input taxa: %d\n - Taxa filtered OUT: %d\n - Taxa RETAINED: %d\n\n", 
              total_initial_taxa, taxa_removed, taxa_kept))
  
  if (nrow(taxa_coords_filtered) == 0) {
    warning("Filtering thresholds are too strict! No taxa passed. Bypassing filter to prevent script crash.")
    taxa_coords_filtered <- taxa_coords
  }
} else {
  message("Taxon loading filter is disabled. Processing all available taxa.")
  taxa_coords_filtered <- taxa_coords
}

# Calculate loading strength and select top N taxa
taxa_coords_filtered$strength_of_taxon_loading <- sqrt(taxa_coords_filtered$Dim1^2 + taxa_coords_filtered$Dim2^2)

top_taxa_arrows <- taxa_coords_filtered %>%
  dplyr::arrange(desc(strength_of_taxon_loading)) %>%
  head(top_asv_n)

# Scale taxa vectors to ~65% max extent
if (nrow(top_taxa_arrows) > 0) {
  scale_factor_taxa    <- (max_point_limit / max(abs(c(top_taxa_arrows$Dim1, top_taxa_arrows$Dim2)))) * 0.65
  top_taxa_arrows$Dim1 <- top_taxa_arrows$Dim1 * scale_factor_taxa
  top_taxa_arrows$Dim2 <- top_taxa_arrows$Dim2 * scale_factor_taxa
}


# =====================================================================
# 4. Final Plot Construction
# =====================================================================
# Base Plot
p_beta_PCA <- ggplot(ordination_df, aes(x = Dim1, y = Dim2, 
                                        color = .data[[color_var]], 
                                        fill  = .data[[color_var]], 
                                        shape = .data[[shape_var]])) + 
  geom_point(size = 4, alpha = 0.8)

# Grouping Ellipses
if (show_ellipses) {
  p_beta_PCA <- p_beta_PCA + 
    stat_ellipse(aes(group = .data[[group_clustering]], colour = .data[[group_clustering]], fill = .data[[group_clustering]]), 
                 geom = "polygon", alpha = 0.1, type = "norm", level = 0.95, linewidth = 0.2)
}

# Environmental Arrows
if (nrow(sig_env_arrows) > 0) {
  p_beta_PCA <- p_beta_PCA +
    geom_segment(data = sig_env_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE, linewidth = 0.6) +
    geom_text_repel(data = sig_env_arrows, aes(x = Dim1, y = Dim2, label = Variable),
                    color = "black", size = 4, fontface = "bold", inherit.aes = FALSE, box.padding = 0.3)
}

# Environmental Centroids
if (nrow(env_centroids) > 0) {
  p_beta_PCA <- p_beta_PCA +
    geom_label_repel(data = env_centroids, 
                     aes(x = Dim1, y = Dim2, label = Clean_Label),
                     color = "black", 
                     fill = "white",
                     size = 3.5, 
                     fontface = "bold.italic", 
                     alpha = 0.9,
                     inherit.aes = FALSE, 
                     box.padding = 0.4,
                     label.padding = 0.2,
                     segment.color = "grey50",
                     segment.size = 0.4)
}

# Top Taxa Arrows
if (nrow(top_taxa_arrows) > 0) {
  p_beta_PCA <- p_beta_PCA +
    geom_segment(data = top_taxa_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "darkblue", alpha = 0.4, inherit.aes = FALSE) +
    geom_text_repel(data = top_taxa_arrows, aes(x = Dim1, y = Dim2, label = Taxon_Label),
                    color = "darkblue", size = 3.5, fontface = "italic", inherit.aes = FALSE,
                    max.overlaps = 20, box.padding = 0.5, point.padding = 0.2)
}

# Optional Sample Labels
if (show_labels) {
  p_beta_PCA <- p_beta_PCA + 
    geom_text_repel(aes(label = .data[[labels]]), size = 3, max.overlaps = 15, show.legend = FALSE)
}

# Aesthetics & Theme
p_beta_PCA <- p_beta_PCA + 
  theme_bw() +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 4, 9, 10, 18, 19, 20)) +
  labs(
    title    = paste(toupper(beta_metric), "Upgraded PCA (Aitchison)"),
    subtitle = paste("PERMANOVA p-value:", format.pval(beta_p_val, digits = 3)),
    x        = paste0("PC1 (", pc1_var, "%)"),
    y        = paste0("PC2 (", pc2_var, "%)"),
    caption  = "Black = Sig. Env Factors | Blue = Top Contributing Taxa"
  ) +
  theme(
    aspect.ratio  = 1, 
    legend.position = "right", 
    text          = element_text(size = 12),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# =====================================================================
# 5. Save Output
# =====================================================================
ggsave(paste0("results/all_wetlands/PCA_Upgraded_", beta_metric, ".png"), 
       p_beta_PCA, width = 12, height = 8, dpi = 300)

message("Upgraded PCA plot generated successfully.")