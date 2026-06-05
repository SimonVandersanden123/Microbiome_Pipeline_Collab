# scripts/05_beta_advanced_vis_PCoA.R
# Based on PCoA_Script_SV foundation
# 1. Extract Coordinates & Merge with Metadata
# -------------------------------------------
# Your PCoA logic: Extract vectors and merge manually with sample_data
ordination_df <- data.frame(ord_beta$vectors[, 1:2]) 
colnames(ordination_df) <- c("Axis.1", "Axis.2")
ordination_df$SampleID <- rownames(ordination_df)
# Extract metadata from current phyloseq object
metadata_current <- data.frame(sample_data(ps_beta_input))
metadata_current$SampleID <- rownames(metadata_current)
# Merge
ordination_df <- merge(ordination_df, metadata_current, by = "SampleID")
# 2. Environmental Fitting (Arrows)
# -------------------------------------------
# Prepare environmental data for significant vector calculation
env_data <- metadata_current %>%
  dplyr::select(dplyr::all_of(env_variables)) %>%
  mutate(across(everything(), function(x) as.numeric(as.character(x))))
         
# Define the 
pcoa_coords <- as.data.frame(ord_beta$vectors[,1:2])
# The function 'envfit'fits environmental vectors or factors onto an ordination. 
# The projections of points onto vectors have maximum correlation with corresponding 
# environmental variables, and the factors show the averages of factor levels.
# For continuous variables this is equal to fitting a linear trend surface (plane in 2D) for a variable (see ordisurf)
# this trend surface can be presented by showing its gradient (direction of steepest increase) using an arrow.
# The environmental variables are the dependent variables that are explained by the ordination scores, 
# and each dependent variable is analysed separately.
# You need a an ordination object or other structure from which the ordination scores can be extracted (including a data frame or matrix of scores).
# Data frame, matrix or vector of environmental variables. The variables can be of mixed type (factors, continuous variables) in data frames.
ef_wunifrac <- vegan::envfit(pcoa_coords, env_data, permutations = 999, na.rm = TRUE)
str(ef_wunifrac)
# Extract and scale arrows
ef_arrows <- as.data.frame(vegan::scores(ef_wunifrac, "vectors")) * vegan::ordiArrowMul(ef_wunifrac)
ef_arrows$Variable <- rownames(ef_arrows)
colnames(ef_arrows)[1:2] <- c("Dim1", "Dim2")

# Filter for significant variables (using your p < 0.05 logic)
sig_arrows <- ef_arrows[ef_wunifrac$vectors$pvals < 0.05, ]
# Manual scaling factors to fix the "too big" arrows
#arrow_reduction_env  <- 0.1, this can now be specified in the rmd file  # Adjust this to shrink black arrows
# Fix Environmental Arrow Scaling as well
if (nrow(sig_arrows) > 0) {
  scale_factor_env <- (max_point_limit / max(abs(c(sig_arrows$Dim1, sig_arrows$Dim2)))) * 0.01
  sig_arrows$Dim1 <- sig_arrows$Dim1 * scale_factor_env
  sig_arrows$Dim2 <- sig_arrows$Dim2 * scale_factor_env
}
# 3. Taxa Correlation (Top ASVs)
# -------------------------------------------
# Unlike PCA (Aitchison), PCoA is based on a distance matrix (e.g., UniFrac) and 
# does not have native species scores (loadings). To see which taxa drive the 
# separation, we use 'envfit' to correlate the original ASV abundances with 
# the fixed PCoA coordinates. This identifies the taxa most strongly 
# associated with the community shifts seen on the axes.

otu_tab <- as(otu_table(ps_beta_input), "matrix")
if(taxa_are_rows(ps_beta_input)) { otu_tab <- t(otu_tab) }

enfit_taxa  <- vegan::envfit(pcoa_coords, otu_tab, permutations = 0)
taxa_scores <- as.data.frame(vegan::scores(enfit_taxa, "vectors"))
taxa_scores$r2 <- enfit_taxa$vectors$r

# --- Get Taxon Names instead of Sequences ---
tax_table_df <- as.data.frame(tax_table(ps_beta_input))
# 1. Simple Merge: Since rows match, just bind the genus column directly
taxa_scores$Taxon_Label <- tax_table_df[[target_level]]
#taxa_scores$Taxon_Label <- tax_table_df$family
# 2. Clean Labels: If Genus is NA, use the ASV rowname instead
taxa_scores$Taxon_Label[is.na(taxa_scores$Taxon_Label)] <- rownames(taxa_scores)[is.na(taxa_scores$Taxon_Label)]
# 3. Scale and Pick Top N

top_taxa_arrows <- taxa_scores %>%
  arrange(desc(r2)) %>%
  head(top_asv_n)

# Auto-scale taxa arrows relative to the data spread
if (nrow(top_taxa_arrows) > 0) {
  scale_factor_taxa <- (max_point_limit / max(abs(c(top_taxa_arrows$Axis.1, top_taxa_arrows$Axis.2)))) * 0.65
  top_taxa_arrows$Dim1 <- top_taxa_arrows$Axis.1 * scale_factor_taxa
  top_taxa_arrows$Dim2 <- top_taxa_arrows$Axis.2 * scale_factor_taxa
}
# 4. Final Plot Construction
# -------------------------------------------
p_beta_PCoA <- ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, 
                                         color = .data[[color_var]], 
                                         fill = .data[[color_var]], 
                                         shape = .data[[shape_var]])) +  
  # Samples Points
  geom_point(size = 4, alpha = 0.8)

# CONDITIONAL: Add Grouping Ellipses
if (show_ellipses) {
  p_beta_PCoA <- p_beta_PCoA + 
    stat_ellipse(aes(group = .data[[group_clustering]]), geom = "polygon", alpha = 0.1, level = 0.95, linewidth = 0.2)
}

# Environmental Arrows (Only Significant)
if (nrow(sig_arrows) > 0) {
  p_beta_PCoA <- p_beta_PCoA +
    geom_segment(data = sig_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE, linewidth = 0.6) +
    geom_text_repel(data = sig_arrows, aes(x = Dim1, y = Dim2, label = Variable),
                    color = "black", size = 4, fontface = "bold", inherit.aes = FALSE, box.padding = 0.3)
}

# Taxa Arrows (Top ASVs)
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

# Styling, Axes Variances, and General Plot Metadata
p_beta_PCoA <- p_beta_PCoA + 
  theme_bw() +
  scale_shape_manual(values = c(7, 9, 15, 0, 16, 1, 17, 2, 18, 5, 19, 6)) + # Safe shapes for over 6 categories
  labs(
    x = paste0("PCoA 1 (", round(ord_beta$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(ord_beta$values$Relative_eig[2] * 100, 1), "%)"),
    title = paste("PCoA Analysis (", toupper(beta_metric), "Distance)"),
    subtitle = paste("PERMANOVA p-value:", format.pval(beta_p_val, digits = 3),
    if(show_labels) paste0("| Points labeled with:", labels) else ""),
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
if(!dir.exists("results")) dir.create("results")
ggsave(paste0("results/PCoA_upgraded_", beta_metric, ".png"), plot = p_beta_PCoA, width = 12, height = 8, dpi = 300)
message("PCoA Upgraded plot successfully generated using SV base logic.")