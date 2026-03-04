# Microbial Data Analysis Pipeline: From Sequences to community traits
**Author:** Simon Vandersanden  
**Contact:** simon.vandersanden@uhasselt.be  
**Last Updated:** 2025-03-24

## 1. Overview
This pipeline transforms raw amplicon sequencing data (ASV/OTU tables) into microbial community traits, which can be inferred to microbial traits or characteristics.We are moving from **counts** of specific bacteria to **functional and phylogenetic traits** that can be used for modelling.

## 2.Environment Setup & Data Containers:
We utilize the phyloseq R package as our primary data container. Think of a phyloseq object as a "relational database" that wraps four distinct tables into one synchronized file. This synchronization is key: if we filter out a "Blank" sample from the metadata, it is automatically removed from the count and taxonomy tables, maintaining data integrity across the pipeline.

### 2.1The Components of the Pipeline
**@otu_table (The Count Matrix)**
This is the heart of the microbial data. In our pipeline, we use ASVs (Amplicon Sequence Variants) generated via the DADA2 algorithm, OTU(Operational Taxonomic Units) could also be used but these are lower resolution. ASVs provide higher resolution by identifying exact biological sequences.
**Rows** (Samples): Represent individual samples (e.g., "nosample1", "MockB").
**Columns** (ASVs): Represent the unique genetic sequences found in the study.
**Values**: Integers representing the "counts" (how many times that specific sequence was detected in that sample).
**Technical Note**: The current data is oriented with taxa_are_rows: FALSE. This means the matrix is structured as:
36 samples X 1346 ASVs
**@tax_table (The identity Table)**
This table acts as the "ID card" for every sequence found in the OTU table. It maps the long DNA strings (ASVs) to a standardized taxonomic hierarchy, allowing us to "name" the bacteria. This information was previously added to the sequences in the DADA2 pipeline during the
**Rows** (Samples): Represent individual samples (e.g., "nosample1", "MockB").
**Columns** (ASVs): Represent the unique genetic sequences found in the study.
**Values**: Integers representing the "counts" (how many times that specific sequence was detected in that sample).
**Technical Note**: The current data is oriented with taxa_are_rows: FALSE. This means the matrix is structured as:
36 samples X 1346 ASVs




```r
# Necessary libraries
library(phyloseq)
library(dplyr)
library(picante)    # For Phylogenetic Diversity
library(microbiome) # For Dominance metrics



