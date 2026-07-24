# scripts/05_beta_diversity_RDA.R
# Prepare the data for the RDA ordination

# Extract the CLR-transformed OTU/ASV table (Species matrix)
# vegan expects samples as rows and species as columns
otu_mat <- as(otu_table(mibi_bmr_clr), "matrix")
if (taxa_are_rows(mibi_bmr_clr)) {
  otu_mat <- t(otu_mat)
}
# Ensure samples are ROWS, taxa are COLUMNS for vegan
# Extract the sample metadata (Environmental matrix)
metadata_raw <- as.data.frame(as(sample_data(mibi_bmr_clr), "data.frame"))
######################################################################################################################
# Different variables of interest: naphthalene_PAH ,PAH_total_16,  PAH_total_10, fraction_C10-C12, fraction_C12-C22, total_C10-C40
# EGV3 ,

# Define your variables of interest
discrete_cols   <- c("Timepoint", "Sample_depth_information", "CW_Number", "CW_Number_MF_number")
continuous_cols <- c("iron2", "phosphorus_total", "concentration_BTEX", "sulfate", 
                     "redoxpot", "temperature", "oxygen", "PAH_total_16")

# naphthalene_PAH,
all_target_cols <- c(discrete_cols, continuous_cols)

# Step 1: Extract ONLY target columns into a fresh table
rda_metadata <- metadata_raw[, all_target_cols]
# Trim out the hidden whitespaces, to ensure if there are empty values that these get filtered out properly
rda_metadata[] <- lapply(rda_metadata, function(x) {
  # Trim hidden whitespace padding
  x <- trimws(as.character(x))
  # If a cell is completely empty, replace it with a true R NA
  x[x == "" | x == "NA"] <- NA
  return(x)
})

# Step 2: Drop rows with NA values in our selected columns 
# RDA cannot compute with missing environmental data
samples_to_keep <- complete.cases(rda_metadata)
rda_metadata_clean <- rda_metadata[samples_to_keep, ]
#which one got dropped? : Short lines to find and print exactly which samples were dropped
dropped_samples <- rownames(rda_metadata)[!samples_to_keep]
cat("Dropped samples due to NA values:", if(length(dropped_samples) > 0) paste(dropped_samples, collapse = ", ") else "None", "\n")

# Step 3: Align your compositional OTU table to match the exact same samples
# (This filters out the matching rows from your CLR transformed matrix)
otu_clr_clean <- otu_mat[samples_to_keep, ]

##############################################################################################################################
# Now that the data is selected, we are going to scale them 
# 1. Convert discrete columns explicitly to Factors (Categorical variables)
for (col in discrete_cols) {
  rda_metadata_clean[[col]] <- as.factor(rda_metadata_clean[[col]])
}

# 2. Convert continuous columns to numeric, catching exactly where text/NAs slip in
for (col in continuous_cols) {
  # Convert to character vector first to check values safely
  char_vals <- as.character(rda_metadata_clean[[col]])
  num_vals  <- as.numeric(char_vals)
  
  # Find if this specific loop introduced any NEW NAs that weren't originally NA
  new_nas_idx <- which(is.na(num_vals) & !is.na(char_vals))
  
  if (length(new_nas_idx) > 0) {
    cat("\n⚠️ Non-numeric text caught in column:", col, "\n")
    # Loop through and print the exact samples and values causing the issue
    for (idx in new_nas_idx) {
      cat("   -> Sample ID:", rownames(rda_metadata_clean)[idx], 
          "| Problematic Value:", paste0("'", char_vals[idx], "'"), "\n")
    }
  }
  
  # Save the converted values back into the table
  rda_metadata_clean[[col]] <- num_vals
}
  
# Scale only the continuous columns (Mean = 0, SD = 1)
rda_metadata_clean[, continuous_cols] <- scale(rda_metadata_clean[, continuous_cols])

# Quick verification check to ensure everything looks correct
str(rda_metadata_clean)
# Plot the correlation matrix to check for which metadata needs to be included in the model
library(corrplot)

cor_matrix <- cor(rda_metadata_clean[, continuous_cols], method = "pearson")
corrplot(cor_matrix, 
         method = "ellipse",       # Represents strength/direction using ellipses
         type = "upper",           # Display only the upper triangle to avoid redundancy
         order = "hclust",         # Hierarchically clusters correlated variables together
         addCoef.col = "black",    # Overlay correlation coefficients as text
         tl.col = "black",         # Text color for labels
         tl.srt = 45,              # Rotate labels 45 degrees
         diag = FALSE)             # Omit the diagonal line (correlation of 1 with self)

# now we split up the metadata table (rda_metadata_clean, into the two seperate matrices)
######################################################################################################
# the modelling starts a global model including all potentiall explanatory variables has to be run as a safety mechanism:
# 1. Define the full global formula
formula_global_loc <- as.formula("otu_clr_clean ~ iron2 + phosphorus_total + concentration_BTEX + 
                                  sulfate + redoxpot + temperature + oxygen + PAH_total_16 + 
                                  Condition(CW_Number_MF_number)")

# 2. Run the global model
model_global_loc <- rda(formula_global_loc, data = rda_metadata_clean)
# 3. Define the permutation constraints (respecting the time series within locations)
ctrl_loc <- how(within = Within(type = "series"), 
                plots  = Plots(strata = rda_metadata_clean$CW_Number_MF_number, type = "none"))
# 4. Run the Global Permutation Test
set.seed(42)
global_test_loc <- anova.cca(model_global_loc, permutations = ctrl_loc)
print(global_test_loc)
# IF and ONLY if this global model is significant, then we can procede to foreward selection.Blanchet’s double-stopping criterion (first criterion)
global_R2_adj <- RsquareAdj(model_global_loc)$adj.r.squared
cat("Global Adjusted R2 ceiling:", global_R2_adj, "\n")
#Blanchet’s double-stopping criterion, second stopping criterion. You need to look at the adjusted R2, and make sure
# the that adjusted R squared of the model you make, must not exceed this R squared.
#####################################################################################################################
# Partial RDA, we will peform a Partial RDA, since our samples are not independant observations, 
# they are linked, either to location or accross time, so we can follow the samples over time.
prda_formula_Sample_Location <- as.formula("otu_clr_clean ~ iron2 + phosphorus_total + concentration_BTEX + 
                    sulfate + redoxpot + temperature + oxygen + PAH_total_16 +
                    Condition(CW_Number_MF_number)")

prda_model_controlled_For_Sample_Location <- rda(prda_formula_Sample_Location, data= rda_metadata_clean)

print(prda_model_controlled_For_Sample_Location)
summary(prda_model_controlled_For_Sample_Location)
# check the variance inflation factors, to see where there is a lot of correlation between the different selected 
# canonical variables , if the vifs are higher than 20, remove the largest one, rerun the model and check them again.
# because the varriance within the one variable  that you remove can influence all other vifs.
vif.cca(prda_model_controlled_For_Sample_Location)
RsquareAdj(prda_model_controlled_For_Sample_Location)

# now we can also control for the varriation that is captured over time
prda_formula_Timepoint <- as.formula("otu_clr_clean ~ iron2 + phosphorus_total + concentration_BTEX + 
                    sulfate + redoxpot + temperature + oxygen + PAH_total_16 +
                    Condition(Timepoint)")

prda_model_controlled_For_Timepoint <- rda(prda_formula_Timepoint, data= rda_metadata_clean)

print(prda_model_controlled_For_Timepoint)
summary(prda_model_controlled_For_Timepoint)

vif.cca(prda_model_controlled_For_Sample_Location)
RsquareAdj(prda_model_controlled_For_Sample_Location)

# We controleren voor Tijd, dus we testen de chemie. 
# De herhaalde metingenstructuur zit nog steeds in CW_Number, dus daar moeten we binnen blokkeren.
ctrl_time <- how(plots = Plots(strata = rda_metadata_clean$CW_Number_MF_number, type = "none"))

# Test algehele significantie van het model
set.seed(42)
anova_overall_time <- anova.cca(prda_model_controlled_For_Timepoint, permutations = ctrl_time)
print(anova_overall_time)

# Test welke chemische variabelen significant zijn als tijd is buitengesloten
anova_terms_time <- anova.cca(prda_model_controlled_For_Timepoint, permutations = ctrl_time, by = "term")
print(anova_terms_time)







# 1. Verification Step: Check VIF again (They should all be < 5 now!)
print("--- New VIF Values ---")
vif.cca(final_bmr_rda)

# 2. Extract Adjusted R-squared (The true variance explained by this clean model)
print("--- Adjusted R-squared ---")
RsquareAdj(final_bmr_rda)

# 3. Test global significance of your new clean metadata model
print("--- Model Significance ---")
anova(final_bmr_rda, permutations = 999)

# Test the unique significance of each term in the model
anova(final_bmr_rda, by = "margin", permutations = 999)


# Create a professional, clean base plot
# scaling = 2 focuses on angles/correlations between variables and samples
plot(final_bmr_rda, scaling = 2, type = "none", 
     main = "Robust Aitchison RDA (Adj. R² = 12.4%, p = 0.001)")

# 1. Plot the samples as small gray dots (keeps the background clean)
points(final_bmr_rda, display = "sites", scaling = 2, pch = 21, bg = "gray80", col = "gray50", cex = 0.8)

# 2. Plot the environmental vector arrows (continuous variables)
text(final_bmr_rda, display = "bp", scaling = 2, col = "darkred", font = 2, cex = 0.9)

# 3. Plot centroids for categorical factors (like depth or CW sections)
text(final_bmr_rda, display = "centroids", scaling = 2, col = "darkblue", font = 2, cex = 0.9)




# Set scaling consistently across all extractions
library(vegan)
library(ggplot2)
library(ggrepel)

sc <- 2

# 1. Extract Sample (Site) Coordinates
sample_scores <- as.data.frame(scores(final_bmr_rda, display = "sites", scaling = sc))
sample_coords <- cbind(sample_scores, rda_metadata_clean)

# FIX: Convert row names (which hold your sample identifiers) into an actual column
sample_coords$Sample_ID <- rownames(sample_coords)

# 2. Extract Continuous Environmental Vectors (Arrows)
all_vectors <- as.data.frame(scores(final_bmr_rda, display = "bp", scaling = sc))
arrow_coords <- all_vectors[rownames(all_vectors) %in% continuous_cols, ]
arrow_coords$Variable <- rownames(arrow_coords)

# 3. Extract Categorical Factor Centroids using "cn"
centroid_scores <- as.data.frame(scores(final_bmr_rda, display = "cn", scaling = sc))
centroid_coords <- centroid_scores
centroid_coords$Factor_Level <- rownames(centroid_coords)
centroid_coords$Factor_Level <- gsub("CW_Number|Sample_depth_information", "", centroid_coords$Factor_Level)

# Visual adjustment factor to stretch arrows/centroids outward for readability
# Visual adjustment factor to stretch arrows/centroids outward for readability
arrow_scale <- 2.5

rda_plot <- ggplot() +
  # A. Origin Crosshairs
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
  
  # B. Microbial Samples (Points) - Colored by Wetland, shaped by Depth
  geom_point(data = sample_coords, 
             aes(x = RDA1, y = RDA2, color = CW_Number, shape = Sample_depth_information), 
             size = 2.8, alpha = 0.75) +
  
  # C. FIXED Layer: Adding Sample Name Labels dynamically using row names
  geom_text_repel(data = sample_coords,
                  aes(x = RDA1, y = RDA2, label = Sample_ID),
                  size = 2.2,                  # Small text so 160 labels fit well
                  color = "gray30",            # Muted gray so it doesn't distract from vectors
                  max.overlaps = 20,           # Drops labels in overly packed spaces to keep clean
                  box.padding = 0.15, 
                  point.padding = 0.1) +
  
  # D. Continuous Chemistry (Arrows)
  geom_segment(data = arrow_coords,
               aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale),
               arrow = arrow(length = unit(0.20, "cm")), 
               color = "firebrick3", linewidth = 0.7, alpha = 0.8) +
  
  # E. Smart Labels for Chemistry Arrows (No overlaps!)
  geom_text_repel(data = arrow_coords,
                  aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = Variable),
                  color = "firebrick4", fontface = "bold", size = 3.8,
                  box.padding = 0.3, force = 2) +
  
  # F. Categorical Feature Centroids (Labels in white boxes)
  geom_label_repel(data = centroid_coords,
                   aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = Factor_Level),
                   color = "dodgerblue4", fill = "white", fontface = "bold", size = 3.5,
                   box.padding = 0.5, force = 3, alpha = 0.9) +
  
  # G. Theme & Aesthetics
  theme_bw() +
  scale_color_brewer(palette = "Set1") + 
  labs(title = "Robust Aitchison RDA Ordination Biplot",
       subtitle = "Model Adj. R² = 12.4%, Global p = 0.001 (All marginal terms p < 0.05)",
       x = "RDA1 (6.78% variance explained)", 
       y = "RDA2 (4.21% variance explained)",
       color = "Constructed Wetland",
       shape = "Sample Depth (m)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold"),
        legend.position = "right")

# Display the final plot

plot(rda_plot)




