# Microbial Data Analysis Pipeline: From Sequences to community traits
**Author:** Simon Vandersanden  
**Contact:** simon.vandersanden@uhasselt.be  
**Last Updated:** 2025-03-24

##  Overview

This pipeline transforms raw amplicon sequencing data (ASV/OTU tables) into microbial community traits, which can be inferred to microbial traits or characteristics.We are moving from **counts** of specific bacteria to **functional and phylogenetic traits** that can be used for modelling.

## 1.Environment Setup & Data Containers:

We utilize the phyloseq R package as our primary data container. Think of a phyloseq object as a "relational database" that wraps four distinct tables into one synchronized file. This synchronization is key: if we filter out a "Blank" sample from the metadata, it is automatically removed from the count and taxonomy tables, maintaining data integrity across the pipeline.

### 1.1 The Components of the Pipeline

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
#### 2.2 Implementation in R
To inspect these components, you can use these functions:

```r
Run this using the .Rmd file:
## 1:Inspecting the different parts of the phyloseq object
```
## 2. Data Pre-processing: Filtering and normalisation
Before the analysis, the data must be cleaned and normalized to ensure that differences between the samples are biological, not technical artifacts.

### 2.1 Filtering
#### Some filtering normally also has been done beforehand: (atleast that is how I do it)
**Sequencing depth threshold**
Removing samples with insufficient reads.
**Contaminant filtering**
Using tools like decontam, here you do need an accurate blanc included in the samples.
**Taxonomic filtering**
Removing:
Eukaryota
Mitochondria
Chloroplast
#### Rarefaction_curves to check if sequencing depth was sufficient:

**rarefaction_curves**: Rarefaction curves plot the number of unique species (ASVs) discovered against the total number of sequences (reads) per sample to determine if the sequencing depth was sufficient to capture the site's full diversity. For data analysis, a curve that reaches a plateau indicates a reliable "steady-state" representation of the community, whereas a steep, rising curve suggests the diversity is under-sampled and the "true" richness remains unknown. Also samples with very low reads are candidates to remove
#### Filtering that is normal to start the analysis of 'clean ready to go phyloseq objects'

##### Low count or prevalance based filtering is done:

**Removal of singletons**: Remove the features that only occur once in one sample, this will be peformed standard.

We apply filters based on the specific research question. This includes:
**Prevalence Filtering**: Removing ASVs that appear in only one or two samples, as these may be sequencing artifacts.

Example: Keep ASVs present in ≥5% of samples, Or present in ≥3 samples

This will : Reduces sparsity, removes potential artifacts, improves statistical power.

**Abundance Filtering**: Removing the sequences which have a low relative abundance

Example: Remove ASVs with total counts < 10, Remove ASVs with relative abundance < 0.01%

##### Then you can filter based on low variance :
Features that remain nearly constant across the experimental conditions are unlikely to be associated with the conditions being studied. Their variability can be assessed using measures such as the interquartile range (IQR). This is not working yet, need to ask at DSI for input.

##### Sample removal or subset selection
Sample removal is usually based on: **Subset Selection**: This is not strictly a "filter" — it’s more: Stratification, Experimental subsetting.

Example: Compare “Contaminated” vs “Reference”, Analyze only PAH polluted samples

So this step is more analytical selection rather than noise filtering.

#### 2.1 Implementation in R
To peform these filtering steps you can use this script:

```r
Run this using the .Rmd file:
## 2:Data Pre-processing: Filtering
```

### 2.2 Normalization & The "Zero Problem"
Because different samples result in different total "reads" (sequencing depth), we cannot compare raw counts directly. Furthermore, microbiome data is **compositional** (counts are parts of a whole), meaning an increase in one taxon's abundance forcedly decreases the relative abundance of others. This gives us a problem when doing statistics on the data:

**Spurious Correlations**: In compositional datasets, traditional correlation coefficients (like Pearson) are biased. Taxa may appear to be negatively correlated simply because they are competing for a fixed 100% of the "pie chart," not because of a biological interaction.

**Differential Abundance Bias**: Without proper normalization, a taxon might appear "significantly different" between treatments solely because the sequencing depth was higher in one group, not because its absolute concentration changed in the soil.

**Beta Diversity**: Metrics like Bray-Curtis are sensitive to library size. Using CLR (Centered Log-Ratio) transformation allows us to move from the "constrained" space of percentages to an "unconstrained" Euclidean space, enabling the use of standard linear models and PCoA.

#### **Handling Sparsity (The Zero Problem)**
Microbial matrices are "sparse" (contain many zeros). Since mathematical transformations like **CLR** involve logarithms, we cannot process zeros directly. We handle this **after filtering** by:

**Pseudo-counts:** Adding a small value (e.g., 1) to all entries.
**Library size based imputation:**: Imputation-based approaches, where the zero’s are imputed differently according to the library size of the sample. Used by LinDA package
**Imputation:** Estimating zero values based on the probability distribution of the detected sequences, various methods exist such as the Bayesian-Multiplicative Replacement (BMR). Bayesian Multiplicative (GBM) model to impute zero values. This method estimates the probability of a taxon being present but undetected based on the sample's total sequencing depth, ensuring that the internal covariance and ratios of the community remain biologically accurate for downstream modeling.
**Matrix Completion**:Matrix completion (often via Probabilistic Matrix Factorization) treats zeros as missing values. It assumes that the data has a low-rank structure—meaning that a few underlying biological factors (like diet, host health, or environment) explain most of the variation in the microbial community. Best for: Predictive modeling.
**Others**: various other imputation appoaches are possible, however often they are included in the packages, such as ANCOMBC2 which uses an other even more advanced type of zero imputation.

**Which imputation method to use depends on the downstream application:**

Use Bayesian-Multiplicative Replacement (BMR) if you are performing standard statistical comparisons (e.g., "Is Taxon A higher in the sick group than the healthy group?"). It is the "gold standard" for preprocessing data before Log-Ratio transformations, which are required to account for the compositional nature of sequencing data.

Use Matrix Completion if you are building predictive machine learning models or if your data is extremely sparse and you believe there is a hidden biological structure (latent factors) that can explain the zeros.
 
#### **Transformation Methods**
**Relative Abundance (TSS):** Normalizes counts to a scale of 0–1 (or 0–100%). While intuitive, it does not solve the compositional bias.
  
**Centered Log-Ratio (CLR):** Transforms the data into a real-number space (Aitchisonian geometry). This removes the "unit sum" constraint and is the gold standard for **linear modeling** and **correlation analysis** in microbial ecology.

##### *When do we use which normalisation method?
**Relative Abundance (TSS):** This is used for simple visualisations, such as simple bar charts, pie charts, simple comparisons (taxon A is twice as abundant in group A compared to group B (this is not statistically differentially abudant, but for simple statements is suffieces). 'Core microbiome analysis, if you want to see which taxa are present in the samples or groups with at least a relative abundance of >0.1%.

**Centered Log-Ratio (CLR):** TThis should be used in statistical modeling, hypothesis testing, and correlation. Basically when you want to use mathematics to interpret or analyse your data and overcome the "compositional" constraint. Correlation analysis, ordination plots, differential abundance testing. **Beta diversity, Correlation Networks, Differential Abundance, Linear Regression & Ordinary differential equation modeling,**

```r
Run this using the .Rmd file:
## 2:Data Pre-processing: Normalisation
```
## 3. Alfa diversity calculations 
Alpha diversity metrics are a general term for metrics that describe the species richness, evenness, or diversity within a sample. Collectively, these metrics contribute to a comprehensive set of traits characterizing the samples, allowing for the determination of key aspects of the microbial community.
### Alfa diversity metrics
Numerous alfa diversity metrics exist. However many different metrics used to determine similar traits yield highly similar results. Therefore a solid selection of the metrics is essential to provide a clear, non-redundant and complete overview of the microbial community traits. 4 main categories exist: richness, dominance, phylogenetic and information metrics (https://doi.org/10.1038/s41598-024-77864-y). 

#### Richness:
Species richness (Observed):A measure of how many different ASV/OTUs are present in each samples.
High richness often suggests a more resilient community.
#### Dominance:
Different measure on the degree in which a taxon or groups of taxa have monopolized the environment. 
##### Evenness 
Describes how fairly the different abundances are distrubuted in the community. In polluted sites, the evenness is typically low since a few specialised species will outcompe the others.
#### Information
Shannon diveristy:
The Shannon Diversity Index (also known as Shannon-Weaver Index or Shannon Entropy) is a measure that quantifies both the richness (total number of different species present) and evenness (how equally abundant the different species are) of species in a community. Unlike simple richness measures that only count the number of different species, the Shannon Index considers how abundance is distributed among those species.

#### 2.1 Implementation in R
```r
Run this using the .Rmd file:
## 3.1:Calculating alfa diversity metrics
## 3.2:To visualise the alfa diversity metric results 
```











