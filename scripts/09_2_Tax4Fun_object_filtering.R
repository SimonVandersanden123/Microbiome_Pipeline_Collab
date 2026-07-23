# Tax4Fun object generation "09_2_Tax4Fun_object_filtering.R"
#load dependencies
library(Tax4Fun2)
# ==========================================
# 1. Load in the data for filtering
# ==========================================
# 1. Read the files into R (skipping the first comment row if Tax4Fun2 includes one)
# Note: check if row 1 is a header or comment; usually Tax4Fun2 uses standard tab tables.
pathway_data <- read.delim(pathway_file, header = TRUE, row.names = 1, check.names = FALSE)
ko_data      <- read.delim(ko_file, header = TRUE, row.names = 1, check.names = FALSE)
# Transpose the data, inorder to maintain a good structure with rows: samples, collumns= functions
pathway_transposed <- t(pathway_data)
head(pathway_transposed)
# ==========================================
# 2. Part to check which groups can be selected
# ==========================================
# 2.1 Extract the level2 and level3 rows to use as our filtering key, this will be used for filtering
level_3_labels <- as.character(pathway_transposed["level3", ])
names(level_3_labels) <- colnames(pathway_data)
level_2_labels <- as.character(pathway_transposed["level2", ])
names(level_2_labels) <- colnames(pathway_data)

# 2.2 Overview of the different pathway levels you can filter on, here you get an overview of the available options:
unique(level_3_labels)
table(level_3_labels)
unique(level_2_labels)
table(level_2_labels)

# ==========================================
# 3. Conditional Filtering Logic
# ==========================================
if (filter_on_KO_level == 'level_3') {
  
  message("Filtering data based on KO Level 3...")
  keep_columns <- which(level_3_labels %in% wanted_categories_level_3)
} else if (filter_on_KO_level == 'level_2') {
  message("Filtering data based on KO Level 2...")
  keep_columns <- which(level_2_labels %in% wanted_categories_level_2)
} else {
  stop("Invalid selection! 'filter_on_KO_level' must be either 'level_2' or 'level_3'.")
}

# ==========================================
# 4. Apply the Filter to Transposed Data
# ==========================================
pathway_transposed <- t(pathway_data)

# Subset the transposed matrix to only keep the selected columns
pathway_filtered <- pathway_transposed[, keep_columns, drop = FALSE]

# Quick check to ensure it worked
message("Dimensions of the filtered object, here you get an overview of the number of KO's maintained:")

print(dim(pathway_filtered))
message("The generated object is called : pathway_filtered")