# Microbiome Analysis Pipeline: From Sequences to Community Traits

This repository contains a modular pipeline designed to transform microbial amplicon data (ASVs) into ecologically relevant traits for downstream data analysis and modeling 

## 📁 Repository Structure

* **`Primary_Documentation.md`**: The primary documentation. Refer to this file for the **scientific theory**, detailed explations behind the different steps of the pipeline, with theoretical reasoning and considerations behind the different steps.
* **`Pipeline_Overview.Rmd`**: The executive summary of the analysis. This file contains a title of the respective steps and pulls code from the `scripts/` folder and generates the visualizations. **Run this to see the results, of the different steps.**
* **`scripts/`**: Modular R scripts. Each script handles one specific step (Loading, Filtering, etc.). These are called by the `.Rmd` via the `source()` command.
* **`data/`**: Contains the input `phyloseq` objects (e.g., `.rds` files).
* **`results/`**: Output directory for generated plots and processed `.csv` tables.

---

## 🚀 How to Use This Pipeline

To reproduce the analysis or explore the data:

1.  **Open the Project**: Double-click `Microbiome_Pipeline_Collab.Rproj` to set the correct working directory.
2.  **Theory First**: If you are new to the microbiome side of this collaboration, please read **`Initial_Version.md`** first to understand the data structures (@otu_table, @tax_table, etc.).
3.  **Run the Analysis**: 
    * Open `docs/Pipeline_Overview.Rmd`.
    * You can now block by block run the different steps of the data analysis
    * This will execute the scripts in the `scripts/` folder in sequence and generate a visual report.

---

## 🛠 Modular Scripts Breakdown

The analysis is broken into numbered steps to ensure a logical flow:

1.  **`01_load_phyloseq.R`**: Imports the `.rds` data and validates the phyloseq object structure.
2.  **`02_filtering.R`**: Filters the data, taking into consideration various principles.
3.  **`03_alpha_diversity.R`**: Calculates within-sample traits (Richness, Shannon, Faith's PD, and Berger-Parker Dominance).
4.  **`04_beta_diversity.R`**: Computes between-sample distances matrices (Bray-Curtis, UniFrac) for community shift modeling.

---
