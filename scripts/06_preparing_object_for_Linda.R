# 06_preparing_object_for_Linda.R

# In order te get a readable plot for LinDA, it is better to agglomerate the reads to a specific phylogenetic level
# they you get information on how mayor groups of bacteria respond to the specified model, 
# otherwise raw sequences will be given.
ps_glommed <- phyloseq::tax_glom(ps_work, taxrank = target_level_linda, NArm = TRUE) #In this step all ASVs which have the same taxonomic level will be agglomerated.
# For plotting we need to change the sequences in the OTU table to a specific taxonomic label otherwise raw sequences will be visible on the plot.
otu_tab <- as.data.frame(t(otu_table(ps_glommed)))# extract the OTU table from the agglomerated object
tax_tab <- as.data.frame(tax_table(ps_glommed)) # extract the taxonomy table from the agglomerated object
#Select the specific phylogenetic level you want to have the labels from:
taxa_ids <- tax_tab[, target_level_linda] #Here the phylogenetic level collumn is extracted, This level of interest
#is then added to the otu and tax table in the following line: 
rownames(otu_tab) <- taxa_ids
rownames(tax_tab) <- taxa_ids
#extract the metadata
meta_df <- data.frame(sample_data(ps_glommed))