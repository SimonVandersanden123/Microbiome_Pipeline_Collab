# scripts/02_rarefaction_curves.R
# Requires: mibi_phylo_object

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is needed for rarefaction curves. Please install it.")
}

message("Calculating rarefaction curves... this may take a moment.")

# Extract OTU table
otu_rare <- as(otu_table(mibi_phylo_object), "matrix")

# Ensure samples are rows for vegan
if (taxa_are_rows(mibi_phylo_object)) {
  otu_rare <- t(otu_rare)
}

# Plot rarefaction curves
vegan::rarecurve(
  otu_rare,
  step = 100,
  cex = 0.5, # Slightly smaller text for cleaner look
  label = FALSE, # Set to FALSE if you have 50+ samples to avoid clutter
  col = "blue",
  main = "Rarefaction Curves"
)