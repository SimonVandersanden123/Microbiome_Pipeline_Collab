library(vegan)

# Extract OTU table
otu <- as(otu_table(mibi_phylo_object), "matrix")

# Ensure samples are rows
if (taxa_are_rows(mibi_phylo_object)) {
  otu <- t(otu)
}

# Plot rarefaction curves
rarecurve(
  otu,
  step = 100,
  cex = 0.6,
  label = TRUE
)
