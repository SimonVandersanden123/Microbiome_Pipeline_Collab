# scripts/05_beta_advanced_vis_PCA.R

# 1. Coordinate & Arrow Extraction (PCA/RDA Style)
# -------------------------------
samp_coords <- as.data.frame(vegan::scores(ord_beta, display = "sites", choices = c(1,2)))
taxa_coords <- as.data.frame(vegan::scores(ord_beta, display = "species", choices = c(1,2)))

colnames(samp_coords) <- c("Dim1", "Dim2")
colnames(taxa_coords) <- c("Dim1", "Dim2")

# Merge sample coords with metadata for plotting
metadata_df <- as(sample_data(ps_beta_input), "data.frame")
samp_coords$SampleID <- rownames(samp_coords)
metadata_df$SampleID <- rownames(metadata_df)
ordination_df <- merge(samp_coords, metadata_df, by = "SampleID")

# 2. Environmental Fitting (Significant Arrows Only)
# -------------------------------
env_data <- metadata_df %>%
  dplyr::select(dplyr::all_of(env_variables)) %>%
  mutate(across(everything(), function(x) as.numeric(as.character(x))))

# PCA objects work directly in envfit
enfit_env <- vegan::envfit(ord_beta, env_data, permutations = 999, na.rm = TRUE)

# Extract, Filter (p < 0.05), and Scale (0.1 reduction like your PCoA)
#arrow_reduction_env <- 0.9 
env_coords <- as.data.frame(vegan::scores(enfit_env, "vectors")) * vegan::ordiArrowMul(enfit_env) * arrow_reduction_env
env_coords$Variable <- rownames(env_coords)
colnames(env_coords)[1:2] <- c("Dim1", "Dim2")
sig_env_arrows <- env_coords[enfit_env$vectors$pvals < 0.05, ]

# 3. Taxa Logic (Top ASVs with Family Names)
# -------------------------------
tax_table_df <- as.data.frame(tax_table(ps_beta_input))
# Scale down taxa arrows (0.3 reduction like your PCoA)
#arrow_reduction_taxa <- 3
#target_level, this is defined in the .RMD file
taxa_coords$Taxon_Label <- tax_table_df[[target_level]]

taxa_coords$Taxon_Label[is.na(taxa_coords$Taxon_Label)] <- rownames(taxa_coords)[is.na(taxa_coords$Taxon_Label)]
taxa_coords$r2 <- sqrt(taxa_coords$Dim1^2 + taxa_coords$Dim2^2) # vector length for PCA

top_taxa_arrows <- taxa_coords %>%
  arrange(desc(r2)) %>%
  head(top_asv_n) %>%
  mutate(Dim1 = Dim1 * arrow_reduction_taxa,
         Dim2 = Dim2 * arrow_reduction_taxa)

# 4. Final Plot Construction
# -------------------------------
p_beta_PCA <- ggplot(ordination_df, aes(x = Dim1, y = Dim2, 
                                        color = .data[[color_var]], 
                                        fill = .data[[color_var]], 
                                        shape = .data[[shape_var]])) + 
  # Grouping Ellipses
  stat_ellipse(aes(group = .data[[group_clustering]]), geom = "polygon", alpha = 0.2, level = 0.95) +
  
  # Environmental Arrows (Black)
  geom_segment(data = sig_env_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = sig_env_arrows, aes(x = Dim1, y = Dim2, label = Variable),
                  color = "black", size = 4, fontface = "bold", inherit.aes = FALSE) +
  
  # Taxa Arrows (Blue)
  geom_segment(data = top_taxa_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue", alpha = 0.4, inherit.aes = FALSE) +
  geom_text_repel(data = top_taxa_arrows, aes(x = Dim1, y = Dim2, label = Taxon_Label),
                  color = "blue", size = 3, fontface = "italic", inherit.aes = FALSE,
                  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.2) +
  
  # Samples
  geom_point(size = 4, alpha = 0.8) + 
  geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 15, show.legend = FALSE) + 
  
  theme_minimal() +
  labs(
    title = paste(toupper(beta_metric), "Pimped PCA (Aitchison)"),
    subtitle = paste("PERMANOVA p-value:", format.pval(beta_p_val, digits = 3)),
    caption = "Black = Sig. Env Factors | Blue = Top Contributing Taxa"
  ) +
  theme(aspect.ratio = 1, legend.position = "right", text = element_text(size = 12))

# 5. Save Output
ggsave(paste0("results/PCA_pimped_", beta_metric, ".png"), p_beta_PCA, width = 10, height = 8)
message("Pimped PCA plot generated.")