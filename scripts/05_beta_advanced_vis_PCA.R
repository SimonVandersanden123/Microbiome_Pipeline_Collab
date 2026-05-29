# scripts/05_beta_advanced_vis_PCA.R

# 1. Coordinate & Arrow Extraction (PCA/RDA Style)
# -------------------------------
samp_coords <- as.data.frame(vegan::scores(ord_beta, display = "sites", choices = c(1,2)))
taxa_coords <- as.data.frame(vegan::scores(ord_beta, display = "species", choices = c(1,2)))

colnames(samp_coords) <- c("Dim1", "Dim2")
colnames(taxa_coords) <- c("Dim1", "Dim2")
# Calculate Variance Explained for Axis Labels
pca_summary <- summary(eigenvals(ord_beta))
pc1_var <- round(pca_summary[2, 1] * 100, 1) # Proportion Explained for PC1
pc2_var <- round(pca_summary[2, 2] * 100, 1) # Proportion Explained for PC2


# Merge sample coords with metadata for plotting
metadata_df <- as(sample_data(ps_beta_input), "data.frame")
samp_coords$SampleID <- rownames(samp_coords)
metadata_df$SampleID <- rownames(metadata_df)
ordination_df <- merge(samp_coords, metadata_df, by = "SampleID")
# Determine max limits of points to automatically auto-scale arrows nicely
max_point_limit <- max(abs(c(ordination_df$Dim1, ordination_df$Dim2)))

# 2. Environmental Fitting (Significant Arrows Only)
# -------------------------------
env_data <- metadata_df %>%
  dplyr::select(dplyr::all_of(numeric_env_variables)) %>%
  mutate(across(everything(), function(x) as.numeric(as.character(x))))

# PCA objects work directly in envfit
enfit_env <- vegan::envfit(ord_beta, env_data, permutations = 999, na.rm = TRUE)

# Auto-scale environmental arrows relative to the data spread
env_coords <- as.data.frame(vegan::scores(enfit_env, "vectors")) 
env_coords$Variable <- rownames(env_coords)
colnames(env_coords)[1:2] <- c("Dim1", "Dim2")
# Filter out the significant variable names directly from the Permanova
permanova_df <- as.data.frame(permanova_marginal) %>% 
  tibble::rownames_to_column("Variable")
# Identify variables where p-value is strictly less than 0.05
sig_permanova_vars <- permanova_df %>%
  filter(`Pr(>F)` < 0.01) %>%
  pull(Variable)
# Track dropped variables for the reporting console output
dropped_vars <- setdiff(permanova_df$Variable, sig_permanova_vars)
# Filter your plot coordinates to ONLY include these globally significant variables
sig_env_arrows <- env_coords[env_coords$Variable %in% sig_permanova_vars, ]

cat(sprintf(
  "\n[Environmental Vector Filter Report (p < 0.01)]:\n - Total input variables : %d\n - Vectors RETAINED     : %d (%s)\n - Vectors DROPPED      : %d (%s)\n\n",
  nrow(permanova_df),
  length(sig_permanova_vars),
  paste(sig_permanova_vars, collapse = ", "),
  length(dropped_vars),
  if(length(dropped_vars) > 0) paste(dropped_vars, collapse = ", ") else "None"
))

# =====================================================================
# Auto-scale environmental arrows relative to the data spread
if (nrow(sig_env_arrows) > 0) {
  # Dynamically scales arrows to occupy roughly 75% of the coordinate space
  scale_factor_env <- (max_point_limit / max(abs(c(sig_env_arrows$Dim1, sig_env_arrows$Dim2)))) * 0.75
  sig_env_arrows$Dim1 <- sig_env_arrows$Dim1 * scale_factor_env
  sig_env_arrows$Dim2 <- sig_env_arrows$Dim2 * scale_factor_env
}
# 3. Taxa Logic (Top ASVs with Family Names)
# -------------------------------
tax_table_df <- as.data.frame(tax_table(ps_beta_input))
# Scale down taxa arrows (0.3 reduction like your PCoA)
taxa_coords$Taxon_Label <- tax_table_df[[target_level]]
taxa_coords$Taxon_Label[is.na(taxa_coords$Taxon_Label)] <- rownames(taxa_coords)[is.na(taxa_coords$Taxon_Label)]
taxa_coords$r2 <- sqrt(taxa_coords$Dim1^2 + taxa_coords$Dim2^2)
total_initial_taxa <- nrow(taxa_coords)
top_taxa_arrows <- taxa_coords %>%
  arrange(desc(r2)) %>%
  head(top_asv_n)


# --- OPTIONAL FILTERING LAYER based on prevalence and abundance---
if (filter_taxon_loadings) {
  message("Applying abundance and prevalence filtering to taxa arrows...")
  ps_for_filter <- mibi_tss # extract the object which is already normalised using TSS
  otu_filter_mat <- as(phyloseq::otu_table(ps_for_filter), "matrix")
  if (!phyloseq::taxa_are_rows(ps_for_filter)) { otu_filter_mat <- t(otu_filter_mat) }
  # 3. Calculate True Sample Metrics
  asv_means      <- rowMeans(otu_filter_mat)
  asv_prevalences <- rowSums(otu_filter_mat > 0) / ncol(otu_filter_mat)
  # 4. Filter based on your custom YAML configurations
  keep_asvs <- names(asv_means[asv_means >= filter_taxon_abundance_filter & 
                                 asv_prevalences >= filter_taxon_prevalence_filter])
  # 5. Subset your taxa coordinates dataframe
  taxa_coords_filtered <- taxa_coords[rownames(taxa_coords) %in% keep_asvs, ]
  # Print an overview of the filtering
  taxa_kept <- nrow(taxa_coords_filtered)
  taxa_removed <- total_initial_taxa - taxa_kept
  cat(sprintf("\n[Taxa Filter Report]:\n - Total input taxa: %d\n - Taxa filtered OUT: %d\n - Taxa RETAINED: %d\n\n", 
              total_initial_taxa, taxa_removed, taxa_kept))
  
  # Safety net: If filtering is too aggressive and leaves 0 taxa, throw a warning and bypass
  if (nrow(taxa_coords_filtered) == 0) {
    warning("Filtering thresholds are too strict! No taxa passed. Bypassing filter to prevent script crash.")
    taxa_coords_filtered <- taxa_coords
  }
} else {
  message("Taxon loading filter is disabled (false). Processing all available taxa.")
  taxa_coords_filtered <- taxa_coords
 
}
# --- CALCULATE LENGTHS & EXTRACT TOP ASVS ---
taxa_coords_filtered$r2 <- sqrt(taxa_coords_filtered$Dim1^2 + taxa_coords_filtered$Dim2^2)
top_taxa_arrows <- taxa_coords_filtered %>%
  arrange(desc(r2)) %>%
  head(top_asv_n)

if (nrow(top_taxa_arrows) > 0) {
  # Auto-scale taxa arrows relative to the data spread
  scale_factor_taxa <- (max_point_limit / max(abs(c(top_taxa_arrows$Dim1, top_taxa_arrows$Dim2)))) * 0.65
  top_taxa_arrows$Dim1 <- top_taxa_arrows$Dim1 * scale_factor_taxa
  top_taxa_arrows$Dim2 <- top_taxa_arrows$Dim2 * scale_factor_taxa
}


# 4. Final Plot Construction
# -------------------------------
# 1. Base Plot Creation (Aesthetics & Points)
p_beta_PCA <- ggplot(ordination_df, aes(x = Dim1, y = Dim2, 
                                        color = .data[[color_var]], 
                                        fill = .data[[color_var]], 
                                        shape = .data[[shape_var]])) + 
  # Samples Points
  geom_point(size = 4, alpha = 0.8)

# 2. CONDITIONAL: Add Grouping Ellipses
if (show_ellipses) {
  p_beta_PCA <- p_beta_PCA + 
    stat_ellipse(aes(group = .data[[group_clustering]]), geom = "polygon", alpha = 0.1, level = 0.95, linewidth = 0.2)
}

# 3. Environmental Arrows (Black) - Unchanged
if (nrow(sig_env_arrows) > 0) {
  p_beta_PCA <- p_beta_PCA +
    geom_segment(data = sig_env_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE, linewidth = 0.6) +
    geom_text_repel(data = sig_env_arrows, aes(x = Dim1, y = Dim2, label = Variable),
                    color = "black", size = 4, fontface = "bold", inherit.aes = FALSE, box.padding = 0.3)
}

# 4. Taxa Arrows (Blue) - Unchanged
if (nrow(top_taxa_arrows) > 0) {
  p_beta_PCA <- p_beta_PCA +
    geom_segment(data = top_taxa_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "darkblue", alpha = 0.4, inherit.aes = FALSE) +
    geom_text_repel(data = top_taxa_arrows, aes(x = Dim1, y = Dim2, label = Taxon_Label),
                    color = "darkblue", size = 3.5, fontface = "italic", inherit.aes = FALSE,
                    max.overlaps = 20, box.padding = 0.5, point.padding = 0.2)
}

# 5. CONDITIONAL: Sample Labeling Logic
if (show_labels) {
  p_beta_PCA <- p_beta_PCA + 
    geom_text_repel(aes(label = .data[[labels]]), size = 3, max.overlaps = 15, show.legend = FALSE)
}

# 6. Styling and General Plot Metadata
# Styling, Axes Variances, and General Plot Metadata
p_beta_PCA <- p_beta_PCA + 
  theme_bw() + # Cleaner grid presentation than minimal for multi-axis plots
  # Manually supply enough robust, easy-to-read shapes (R shape codes 15 to 25 are filled/clear)
  scale_shape_manual(
    values = c(16, 17, 15, 3, 7, 8, 4, 9, 10, 18, 19, 20)
  ) +
  labs(
    title = paste(toupper(beta_metric), "Upgraded PCA (Aitchison)"),
    subtitle = paste("PERMANOVA p-value:", format.pval(beta_p_val, digits = 3)),
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    caption = "Black = Sig. Env Factors | Blue = Top Contributing Taxa"
  ) +
  theme(
    aspect.ratio = 1, 
    legend.position = "right", 
    text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# 5. Save Output
ggsave(paste0("results/PCA_Upgraded_", beta_metric, ".png"), p_beta_PCA, width = 12, height = 8, dpi = 300)
message("Upgraded PCA plot generated.")