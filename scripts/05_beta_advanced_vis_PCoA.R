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
sig_arrows <- sig_arrows %>%
  mutate(Dim1 = Dim1 * arrow_reduction_env,
         Dim2 = Dim2 * arrow_reduction_env)
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
# Note: Added 0.3 scaling here to stop them from being "too big"
arrow_scale <- vegan::ordiArrowMul(enfit_taxa) * arrow_reduction_env

top_taxa_arrows <- taxa_scores %>%
  arrange(desc(r2)) %>%
  head(top_asv_n) %>%
  mutate(Dim1 = Axis.1 * arrow_scale,
         Dim2 = Axis.2 * arrow_scale)

# 4. Final Plot Construction
# -------------------------------------------
p_beta_PCoA <- ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, 
                                         color = .data[[color_var]], 
                                         fill = .data[[color_var]], 
                                         shape = .data[[shape_var]])) +  
  # Grouping Ellipses
  stat_ellipse(aes(group = .data[[group_clustering]]), geom = "polygon", alpha = 0.2, level = 0.95) + 
  
  # Environmental Arrows (Only Significant)
  geom_segment(data = sig_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = sig_arrows, aes(x = Dim1, y = Dim2, label = Variable),
                  color = "black", size = 4, fontface = "bold", inherit.aes = FALSE) +
  
  # Taxa Arrows (Top ASVs)
  geom_segment(data = top_taxa_arrows, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue", alpha = 0.4, inherit.aes = FALSE) +
  geom_text_repel(data = top_taxa_arrows, aes(x = Dim1, y = Dim2, label = Taxon_Label),
                  color = "blue", size = 3, fontface = "italic", inherit.aes = FALSE,max.overlaps = Inf,    # Forces all labels to show
                  box.padding = 0.5,     # Gives labels more room to move
                  point.padding = 0.2) + # Distance from the arrow tip 
  
  # Samples
  geom_point(size = 4, alpha = 0.8) + 
  geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 15, show.legend = FALSE) +  
  
  theme_minimal() +
  labs(x = paste0("PCoA 1 (", round(ord_beta$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA 2 (", round(ord_beta$values$Relative_eig[2] * 100, 2), "%)"),
       title = paste("PCoA Analysis (", beta_metric, "Distance)"),
       subtitle = "Black = Sig. Env Factors | Blue = Top Contributing Taxa") +
  theme(legend.position = "right", text = element_text(size = 12), aspect.ratio = 1)

# 5. Save Output
if(!dir.exists("results")) dir.create("results")
ggsave(paste0("results/PCoA_pimped_", beta_metric, ".png"), plot = p_beta_PCoA, width = 10, height = 8)

message("PCoA Pimped plot successfully generated using SV base logic.")