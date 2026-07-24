library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
tax4fun_pathway_data_selection_matrix_bmr_clr  <-readRDS("results/all_wetlands/tax4fun_pathway_data_selection_matrix_bmr_clr.rds")

# now we have the CLR transformed data, so we can now do a PCA, so an unconstrained RDA
ord_tax4fun_traits <- vegan::rda(tax4fun_pathway_data_selection_matrix_bmr_clr)
biplot(ord_tax4fun_traits)
str(tax4fun_pathway_data)
# 1. Extract the exact variance explained for the axes
eigenvalues <- summary(ord_tax4fun_traits)$cont$importance
pc1_var <- round(eigenvalues[2, 1] * 100, 1) 
pc2_var <- round(eigenvalues[2, 2] * 100, 1) 

# 2. Extract Sample (Sites) coordinates
sample_scores <- as.data.frame(scores(ord_tax4fun_traits, display = "sites"))

# FIX: Map the TRUE text SampleID column directly from your source data, matching the row order
sample_scores$SampleID <- tax4fun_pathway_data$SampleID

# 3. MERGE your metadata into the PCA coordinates (This will now match perfectly!)
metadata_sub <- tax4fun_pathway_data[, c("SampleID", "short_name", "wetland_sample", "CW_number", "CW_section", "DNA_concentration", "timepoint", "sample_depth_information","MF_number")]
sample_scores <- left_join(sample_scores, metadata_sub, by = "SampleID")

# 4. Extract Vector (Species/Traits) coordinates
vector_scores <- as.data.frame(scores(ord_tax4fun_traits, display = "species"))
vector_scores$Pathway <- rownames(vector_scores)

# 5. Build the customized colored PCA
ggplot() +
  # Draw the sample points, colored by your depth metadata
  geom_point(data = sample_scores, aes(x = PC1, y = PC2, color = sample_depth_information), size = 3, alpha = 0.8) +
  
  # Add sample text labels using short_name
  geom_text_repel(data = sample_scores, aes(x = PC1, y = PC2, label = short_name),
                  size = 2, color = "grey40", max.overlaps = 20) + # Increased max.overlaps to stop warnings
  
  # Draw the pathway vector arrows
  geom_segment(data = vector_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")), color = "darkred", alpha = 0.3, linewidth = 0.4) +
  
  # Add labels to the vector arrows
  geom_text_repel(data = vector_scores, aes(x = PC1, y = PC2, label = Pathway),
                  color = "darkred", size = 2.5, fontface = "bold", max.overlaps = 20) +
  
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "Functional PCA Colored by Sample Depth",
    color = "Sample Depth"
  ) +
  
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )