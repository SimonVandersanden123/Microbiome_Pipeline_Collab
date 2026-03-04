# Microbial Data Analysis Pipeline: From Sequences to community traits
**Author:** Simon Vandersanden  
**Contact:** simon.vandersanden@uhasselt.be  
**Last Updated:** 2025-03-24

## 1. Overview

This pipeline transforms raw amplicon sequencing data (ASV/OTU tables) into microbial community traits, which can be inferred to microbial traits or characteristics.We are moving from **counts** of specific bacteria to **functional and phylogenetic traits** that can be used for modelling.

## 2.Environment Setup & Data Containers:

We utilize the phyloseq R package as our primary data container. Think of a phyloseq object as a "relational database" that wraps four distinct tables into one synchronized file. This synchronization is key: if we filter out a "Blank" sample from the metadata, it is automatically removed from the count and taxonomy tables, maintaining data integrity across the pipeline.

### 2.1 The Components of the Pipeline

#### **@otu_table (The Count Matrix)**

This is the heart of the microbial data. In our pipeline, we use ASVs (Amplicon Sequence Variants) generated via the DADA2 algorithm, OTU(Operational Taxonomic Units) could also be used but these are lower resolution. ASVs provide higher resolution by identifying exact biological sequences.

**Rows** (Samples): Represent individual samples (e.g., "nosample1", "MockB").

**Columns** (ASVs): Represent the unique genetic sequences found in the study.

**Values**: Integers representing the "counts" (how many times that specific sequence was detected in that sample).

**Technical Note**: The current data is oriented with taxa_are_rows: FALSE. This means the matrix is structured as:
36 samples X 1346 ASVs

#### **@tax_table (The identity Table)**

This table acts as the "ID card" for every sequence found in the OTU table. It maps the long DNA strings (ASVs) to a standardized taxonomic hierarchy, allowing us to "name" the bacteria.

**Structure & Taxonomic Levels**:
The table is a matrix where each row corresponds to a unique ASV, and the columns represent the following hierarchical ranks (from broadest to most specific): **domain**,**phylum**,**class**,**order**,**family**,**genus**,**species**, 

**Example Hierarchy (Human)**: Eukarya → Chordata → Mammalia → Primates → Hominidae → Homo sapiens

**Modeling Utility**:

**1.Complexity Reduction**: Microbial datasets often contain thousands of unique ASVs. By aggregating counts at the Family or Order level, we reduce the dimensionality of the data, making it more manageable for statistical models without losing significant biological signals.

**2. Trait Inference**: Certain taxonomic groups are strongly associated with specific metabolic traits. For example, finding a high abundance of Sphingomonadaceae often indicates a high potential for PAH degradation.

**3. Functional Prediction**:This table is the primary input for tools like Tax4Fun or PICRUSt2. These packages use the taxonomic identity to "link" our sequences to genomic databases, allowing us to predict the functional gene profile (e.g., presence of dioxygenase genes) of the community without performing expensive shotgun metagenomics.

#### **@sam_data (The Metadata/Context)**

This is where the microbiology meets the chemistry. It contains all environmental variables related to the samples, but also arbitrary grouping you can impose on the groups, like locations, or timepoints,...

**typical variables**: PAH concentration, level of electron acceptors, sample groups, timepoints,....

**format** a standard data frame where the rows correspond exactly to the Sample IDs found in the @otu_table.

#### **@phy_tree (The Phylogenetic Tree)**

A phylogenetic tree represents the "genetic distance" between the different sequences. TLDR: sequences that look alike, are likely more related, so they are closer in the phylogenetic tree. 

**Use**: This is often optional but can be used for calculating more advanced microbial community metrics such as **Faith's Phylogenetic Diversity** or **UniFrac distances** which take into account how closely related the different species are. So it is a metric that takes into account how similar the sequences are, instead of just looking at 'how many' different sequences that are found. Because this might be important in the biological question you are asking.
## Implementation in R
To inspect these components, you can use these functions:

```r
# Load in the necessary libraries:
library(phyloseq)
library(dplyr)
library(picante)    # For Phylogenetic Diversity
library(microbiome) # For Dominance metrics
# Look at the data in your Phyloseq object:
# Load in the phyloseq object, which was generated at the end of the DADA2 pipeline.
Phylo_Object <- readRDS("ps_IdTax_noncontam_NoMito_No_Chloro.rds")

# To see the count data (OTU table)
otu_table_data <- as(otu_table(Phylo_Object), "matrix")
head(otu_table_data[, 1:5]) # View first 5 ASVs

# To see the Taxonomy
tax_table_data <- as(tax_table(Phylo_Object), "matrix")
head(tax_table_data)

# To see the Metadata
metadata_data <- as(sample_data(Phylo_Object), "data.frame")
head(metadata_data)









