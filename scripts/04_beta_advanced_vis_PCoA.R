# ==============================================================================
# scripts/05_beta_advanced_vis_PCoA.R
# Refactored to match the upgraded PCA workflow architecture
# ==============================================================================

# 1. Coordinate Extraction & Metadata Merging
# ------------------------------------------------------------------------------
# Extract the first two dimensions from the PCoA ordination object
ordination_df <- as.data.frame(ord_beta$vectors[, 1:2]) 
colnames(ordination_df) <- c("Axis.1", "Axis.2")
ordination_df$SampleID <- rownames(ordination_df)

# Extract and match metadata from the active phyloseq object
metadata_current <- data.frame(sample_data(ps_beta_input))
metadata_current$SampleID <- rownames(metadata_current)

# Merge coordinates with metadata
ordination_df <- merge(ordination_df, metadata_current, by = "SampleID")

# Calculate data spread threshold for dynamic arrow scaling
max_point_limit <- max(abs(c(ordination_df$Axis.1, ordination_df$Axis.2)))


# 2. Environmental Fitting (Vectors & Factors)
# ------------------------------------------------------------------------------
# Extract target numeric and categorical environmental variables
env_data_pre <- metadata_df %>%
  dplyr::select(dplyr::all_of(c(numeric_env_variables, categ_env_variables)))

# Ensure environmental dataset perfectly matches the ordination rows
ord_sample_ids <- rownames(ord_beta$vectors) 
env_data <- env_data_pre[ord_sample_ids, , drop = FALSE]

# Extract raw 2D coordinates for vector fitting
pcoa_coords <- as.data.frame(ord_beta$vectors[, 1:2])

# Fit environmental variables onto the PCoA space
ef <- vegan::envfit(pcoa_coords, env_data, permutations = 999, na.rm = TRUE)

# Safe extraction of vector scores (Guards against 0 numeric vectors matched)
vector_scores <- vegan::scores(ef, "vectors")

if (!is.null(vector_scores) && nrow(vector_scores) > 0) {
  ef_arrows <- as.data.frame(vector_scores) * vegan::ordiArrowMul(ef)
  ef_arrows$Variable <- rownames(ef_arrows)
  colnames(ef_arrows)[1:2] <- c("Dim1", "Dim2")
  
  sig_arrows <- ef_arrows[ef$vectors$pvals < 0.05, , drop = FALSE]
  
  if (nrow(sig_arrows) > 0) {
    scale_factor_env <- (max_point_limit / max(abs(c(sig_arrows$Dim1, sig_arrows$Dim2)))) * 0.75
    sig_arrows$Dim1  <- sig_arrows$Dim1 * scale_factor_env
    sig_arrows$Dim2  <- sig_arrows$Dim2 * scale_factor_env
  }
} else {
  sig_arrows <- data.frame()
}

# 3. Taxa Correlation & Optional Abundance/Prevalence Filtering
# ------------------------------------------------------------------------------
otu_tab <- as(otu_table(ps_beta_input), "matrix")
if (taxa_are_rows(ps_beta_input)) { 
  otu_tab <- t(otu_tab) 
}

enfit_taxa  <- vegan::envfit(pcoa_coords, otu_tab, permutations = 0)
taxa_scores <- as.data.frame(vegan::scores(enfit_taxa, "vectors"))
taxa_scores$r2 <- enfit_taxa$vectors$r
colnames(taxa_scores)[1:2] <- c("Dim1", "Dim2") # Standardize axis column names

tax_table_df <- as.data.frame(tax_table(ps_beta_input))
taxa_scores$Taxon_Label <- tax_table_df[[target_level]]
taxa_scores$Taxon_Label[is.na(taxa_scores$Taxon_Label)] <- rownames(taxa_scores)[is.na(taxa_scores$Taxon_Label)]

total_initial_taxa <- nrow(taxa_scores)

if (isTRUE(filter_taxon_loadings)) {
  message("Applying abundance and prevalence filtering to PCoA taxa vectors...")
  
  ps_for_filter  <- mibi_tss 
  otu_filter_mat <- as(phyloseq::otu_table(ps_for_filter), "matrix")
  if (!phyloseq::taxa_are_rows(ps_for_filter)) { otu_filter_mat <- t(otu_filter_mat) }
  
  asv_means       <- rowMeans(otu_filter_mat)
  asv_prevalences <- rowSums(otu_filter_mat > 0) / ncol(otu_filter_mat)
  
  n_samples                    <- nrow(metadata_complete)
  min_sample_count             <- config$Beta_Diversity$advanced$filter_taxon_prevalence_filter
  calculated_prevalence_filter <- min_sample_count / n_samples
  
  keep_asvs <- names(asv_means[asv_means >= filter_taxon_abundance_filter & 
                                 asv_prevalences >= calculated_prevalence_filter])
  
  taxa_scores_filtered <- taxa_scores[rownames(taxa_scores) %in% keep_asvs, ]
  taxa_kept            <- nrow(taxa_scores_filtered)
  taxa_removed         <- total_initial_taxa - taxa_kept
  
  cat(sprintf("\n[Taxa Filter Report]:\n - Total input taxa: %d\n - Taxa filtered OUT: %d\n - Taxa RETAINED: %d\n\n", 
              total_initial_taxa, taxa_removed, taxa_kept))
  
  if (nrow(taxa_scores_filtered) == 0) {
    warning("Filtering thresholds are too strict! No taxa passed. Bypassing filter to prevent script crash.")
    taxa_scores_filtered <- taxa_scores
  }
} else {
  message("Taxon loading filter is disabled (FALSE). Processing all available taxa.")
  taxa_scores_filtered <- taxa_scores
}

top_taxa_arrows <- taxa_scores_filtered %>%
  dplyr::arrange(desc(r2)) %>%
  head(top_asv_n)

if (nrow(top_taxa_arrows) > 0) {
  scale_factor_taxa    <- (max_point_limit / max(abs(c(top_taxa_arrows$Dim1, top_taxa_arrows$Dim2)))) * 0.75
  top_taxa_arrows$Dim1 <- top_taxa_arrows$Dim1 * scale_factor_taxa
  top_taxa_arrows$Dim2 <- top_taxa_arrows$Dim2 * scale_factor_taxa
}

# 4. Final Plot Construction
# ------------------------------------------------------------------------------
p_beta_PCoA <- ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, 
                                         color = .data[[color_var]], 
                                         fill = .data[[color_var]], 
                                         shape = .data[[shape_var]])) +  
  # Ground sample points layer
  geom_point(size = 4, alpha = 0.8)

# CONDITIONAL: Add Grouping Ellipses
if (show_ellipses) {
  p_beta_PCoA <- p_beta_PCoA + 
    stat_ellipse(aes(group = .data[[group_clustering]]), geom = "polygon", alpha = 0.1, level = 0.95, linewidth = 0.2)
}

# CONDITIONAL: Environmental Arrows (Only Significant)
if (nrow(sig_arrows) > 0) {
  p_beta_PCoA <- p_beta_PCoA +
    geom_segment(data = sig_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE, linewidth = 0.6) +
    geom_text_repel(data = sig_arrows, aes(x = Dim1, y = Dim2, label = Variable),
                    color = "black", size = 4, fontface = "bold", inherit.aes = FALSE, box.padding = 0.3)
}

# CONDITIONAL: Taxa Arrows (Top ASVs)
if (nrow(top_taxa_arrows) > 0) {
  p_beta_PCoA <- p_beta_PCoA +
    geom_segment(data = top_taxa_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue", alpha = 0.4, inherit.aes = FALSE) +
    geom_text_repel(data = top_taxa_arrows, aes(x = Dim1, y = Dim2, label = Taxon_Label),
                    color = "blue", size = 3.5, fontface = "italic", inherit.aes = FALSE,
                    max.overlaps = 20, box.padding = 0.5, point.padding = 0.2)
}

# CONDITIONAL: Sample Labeling Logic
if (show_labels) {
  p_beta_PCoA <- p_beta_PCoA + 
    geom_text_repel(aes(label = .data[[labels]]), size = 3, max.overlaps = 15, show.legend = FALSE)
}

# Theme, Annotations, and Global Aesthetics Layout
p_beta_PCoA <- p_beta_PCoA + 
  theme_bw() +
  scale_shape_manual(values = c(7, 9, 15, 0, 16, 1, 17, 2, 18, 5, 19, 6)) + 
  labs(
    x = paste0("PCoA 1 (", round(ord_beta$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(ord_beta$values$Relative_eig[2] * 100, 1), "%)"),
    title = paste("PCoA Analysis (", toupper(beta_metric), "Distance)"),
    subtitle = paste("PERMANOVA p-value:", format.pval(beta_p_val, digits = 3),
                     if(show_labels) paste0(" | Points labeled with: ", labels) else ""),
    caption = "Black = Sig. Env Factors | Blue = Top Contributing Taxa"
  ) +
  theme(
    aspect.ratio = 1,
    legend.position = "right", 
    text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# 5. Output Management
# ------------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")
ggsave(paste0("results/PCoA_upgraded_", beta_metric, ".png"), plot = p_beta_PCoA, width = fig_width, height = fig_height, dpi = 300)
message("PCoA Upgraded plot successfully generated using SV base logic.")