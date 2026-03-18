# 06_differential_abundance_LinDA.R
# Differential abundance using Linear Models for Microbiome Compositional Data

# 1. Prepare Data from Phyloseq
# -------------------------------------------
# Use the filtered counts (not normalized)
otu_tab <- as.data.frame(otu_table(ps_work))
# LinDA expects features as rows
if(!taxa_are_rows(ps_work)) { otu_tab <- t(otu_tab) }

meta_df <- data.frame(sample_data(ps_work))

# 2. Configuration & Reference Setting
# -------------------------------------------
# This is done in the .Rmd file

# 3. Run LinDA
# -------------------------------------------
# This is done in the .Rmd file
# 4. Extract and Annotate Results
# -------------------------------------------
# Let's look at the interaction effect or the main effect
# Check available variables: print(names(linda_res$output))
var_to_plot <- "Groundwater_PositionBelow_GW" # Example variable

# Create a results table with Taxonomy for the specific variable
res_df <- as.data.frame(linda_res$output[[var_to_plot]])
tax_tab <- as.data.frame(tax_table(ps_work))
res_df <- merge(res_df, tax_tab, by = "row.names")
colnames(res_df)[1] <- "ASV_ID"

# 5. Visualisation
# -------------------------------------------
if(!dir.exists("results/linda")) dir.create("results/linda", recursive = TRUE)

# Generate plots for all variables in the model
plots <- linda.plot(
  linda.obj = linda_res,
  variables.plot = names(linda_res$output),
  titles = names(linda_res$output),
  alpha = 0.05,
  lfc.cut = 1,
  legend = TRUE
)

# Save Volcano Plot for the main variable
ggsave(
  filename = paste0("results/linda/Volcano_", var_to_plot, ".png"),
  plot = plots$plot.volcano[[1]],
  width = 12, height = 9
)

# Save LFC Plot (Log Fold Change)
ggsave(
  filename = paste0("results/linda/LFC_Plot_", var_to_plot, ".png"),
  plot = plots$plot.lfc[[1]],
  width = 10, height = 12
)

message("LinDA analysis complete. Significant taxa found for ", var_to_plot, ": ", sum(res_df$reject))