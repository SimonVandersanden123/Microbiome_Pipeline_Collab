# scripts/05_beta_diversity_RDA.R
# Extract the CLR-transformed OTU/ASV table (Species matrix)
otu_mat <- as(otu_table(mibi_bmr_clr), "matrix")
if (taxa_are_rows(mibi_bmr_clr)) {
  otu_mat <- t(otu_mat)
}
# Ensure samples are ROWS, taxa are COLUMNS for vegan
# Extract the sample metadata (Environmental matrix)
# Ensure samples are ROWS, taxa are COLUMNS for vegan
# Extract the sample metadata (Environmental matrix)
metadata_raw <- as.data.frame(as(sample_data(mibi_bmr_clr), "data.frame"))

# Create your active data frame using your naming convention
metadata_correct_format <- metadata_raw

# Prepare the data for the RDA ordination

# Define your variables of interest
categ_env_variables   <- config$Beta_Diversity$advanced_for_rda$categ_env_variables
numeric_env_variables <- config$Beta_Diversity$advanced_for_rda$numeric_env_variables

# =====================================================================
# 1. AUTOMATED FACTOR CONVERSION LAYER (Generalized)
# =====================================================================
for (cat_var in categ_env_variables) {
  if (cat_var %in% colnames(metadata_correct_format)) {
    # Clean up column to a character format first, then cast to a standard factor
    # This ensures R handles any messy mixed-type inputs or text labels cleanly
    metadata_correct_format[[cat_var]] <- factor(as.character(metadata_correct_format[[cat_var]]))
    message(paste("Successfully converted to factor:", cat_var))
  } else {
    warning(paste("Configured categorical variable not found in dataset columns:", cat_var))
  }
}
# =====================================================================
# 2. AUTOMATED NUMERICAL CONVERSION LAYER (Generalized)
# =====================================================================
for (num_var in numeric_env_variables) {
  if (num_var %in% colnames(metadata_correct_format)) {
    # Step A: Convert to character and strip common hidden spacing/text issues
    clean_chars <- as.character(metadata_correct_format[[num_var]])
    clean_chars <- trimws(clean_chars) # Remove leading/trailing whitespaces
    # Step B: Cast explicitly to numeric double
    # Suppress warnings temporarily so R doesn't flood the console if 
    # text like "missing" or "BDL" naturally turns into NA
    suppressWarnings({
      numeric_vector <- as.numeric(clean_chars)
    })
    # Save the clean numeric vector back into your corrected formatting table
    metadata_correct_format[[num_var]] <- numeric_vector
    message(paste("Successfully converted to numeric metric:", num_var))
    # Optional Safety Check: Report if the conversion introduced any new NAs
    raw_na_count <- sum(is.na(clean_chars)) # check NAs before parsing
    new_na_count <- sum(is.na(numeric_vector))
    if (new_na_count > raw_na_count) {
      warning(paste0(
        "Variable '", num_var, "' contained non-numeric text strings (e.g., 'ND', '<LOD', or blanks) ",
        "which have been automatically forced to NA values."
      ))
    }
  } else {
    warning(paste("Configured numeric variable not found in dataset columns:", num_var))
  }
}
all_target_cols <- c(categ_env_variables, numeric_env_variables)
#-----------------------------------------------------------------------------
# Step 1: Extract ONLY target metadata columns into a fresh table
#-----------------------------------------------------------------------------
# Step 1: Extract ONLY target columns into a fresh table
rda_metadata <- metadata_correct_format[, all_target_cols]

# Step 2: Drop rows with NA values in our selected columns 
# RDA cannot compute with missing environmental data
samples_to_keep <- complete.cases(rda_metadata)
rda_metadata_clean <- rda_metadata[samples_to_keep, ]

# Step 3: Align your compositional OTU table to match the exact same samples
# (This filters out the matching rows from your CLR transformed matrix)
otu_clr_clean <- otu_mat[samples_to_keep, ]

# Scale only the continuous columns (Mean = 0, SD = 1)
rda_metadata_clean[, numeric_env_variables] <- scale(rda_metadata_clean[, numeric_env_variables])

# Quick verification check to ensure everything looks correct
str(rda_metadata_clean)

#-------------------------------------------------------------------------------------
#Check if there are correlations between you different environmental data
#-------------------------------------------------------------------------------------
# Plot the correlation matrix to check for which metadata needs to be included in the model
# 1. Clone your clean numeric data so we don't overwrite your master metadata
cor_data <- rda_metadata_clean[, numeric_env_variables]

# 2. Rename columns to short codes so they don't overlap in the plot
# (Make sure the order matches your numeric_env_variables exactly!)
colnames(cor_data) <- c(
  "BTEX_ug", 
  "PAH_ug", 
  "Total_Fe2_acceptor", 
  "SO4_acceptor", 
  "Mn2_acceptor", 
  "Total_Acceptors_consumed", 
  "Donor_Acceptor_Balance"
)

# 3. Re-calculate the Pearson matrix with the clean, short names
cor_matrix <- cor(cor_data, method = "pearson")

# 4. Set plot margins to give labels breathing room 
# par(mar = c(bottom, left, top, right))
par(mar = c(1, 1, 4, 1)) 

# 5. Run the optimized corrplot
corrplot(cor_matrix, 
         method = "ellipse",       # Ellipses for strength/direction
         type = "upper",           # Upper triangle only
         order = "hclust",         # Group highly correlated variables together
         addCoef.col = "black",    # Text correlation coefficients
         number.cex = 0.8,         # Scale down the number text size slightly so it fits inside ellipses
         tl.col = "black",         # Text color for labels
         tl.srt = 45,              # Rotate labels 45 degrees
         tl.cex = 0.85,            # Scale down label text size so it stays on screen
         diag = FALSE)             # Omit self-diagonal (1.00)

#Foreward model selection:
#
# 1. Define the absolute minimum model (Intercept only - no metadata)
null_model <- rda(otu_clr_clean ~ 1, data = rda_metadata_clean)

# 2. Define the absolute maximum model (All your metadata variables)
full_model <- rda(otu_clr_clean ~ Timepoint + Sample_depth_information + CW_Number + 
                    iron2 + phosphorus_total + concentration_BTEX + sulfate + 
                    redoxpot + temperature + oxygen + PAH_total_16, 
                  data = rda_metadata_clean)

# 3. Run automatic forward selection
# This will step-by-step select variables based on AIC and permutation significance
stepwise_model <- ordistep(null_model, 
                           scope = formula(full_model), 
                           direction = "forward", 
                           permutations = 499)

# 4. Check which variables made the final cut
summary(stepwise_model)
vif.cca(stepwise_model)
# here we saw that temperature and the timepoint are not very good varriables due to collinearity
# e.g. temperature can be predicted for a large part by timepoint.
vif.cca(stepwise_model)

# Here you can specify the model formula to include the different metadata
rda_formula <- as.formula("otu_clr_clean ~ Timepoint + Sample_depth_information + CW_Number + iron2 + phosphorus_total + concentration_BTEX + sulfate + redoxpot + temperature + oxygen + PAH_total_16")

# Pruned formula: Dropping Timepoint to fix temperature collinearity
pruned_formula <- as.formula("otu_clr_clean ~ CW_Number + Sample_depth_information + temperature + 
                             iron2 + phosphorus_total + concentration_BTEX + sulfate + redoxpot + PAH_total_16")
# Re-run the RDA
final_bmr_rda <- rda(pruned_formula, data = rda_metadata_clean)

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




