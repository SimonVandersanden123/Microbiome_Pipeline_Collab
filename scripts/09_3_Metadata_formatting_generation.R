#We are going to stepwisely deconstruct and format the table into a suitable format for copy pasting it into the metadata for the phyloseq 
# Here we aim to extract the information and add it to the phyloseq metadata
# 1. Separate the 54 biological samples from the 3 hierarchical character rows
sample_rows  <- setdiff(rownames(pathway_filtered), c("level1", "level2", "level3"))
level1_names <- as.character(pathway_filtered["level1", ])
level2_names <- as.character(pathway_filtered["level2", ])
level3_names <- as.character(pathway_filtered["level3", ])
ko_numbers   <- colnames(pathway_filtered)
# 2. Extract and convert the text matrix part into a true numeric matrix
raw_numeric_part     <- pathway_filtered[sample_rows, ]
clean_numeric_matrix <- matrix(as.numeric(raw_numeric_part), 
                               nrow = nrow(raw_numeric_part), 
                               ncol = ncol(raw_numeric_part))
rownames(clean_numeric_matrix) <- rownames(raw_numeric_part)
colnames(clean_numeric_matrix) <- colnames(raw_numeric_part)
# 3. Create a clean map framework connecting KO to all its hierarchical levels
hierarchy_map <- data.frame(
  KO_Number = ko_numbers,
  Level1    = level1_names,
  Level2    = level2_names,
  Level3    = level3_names,
  stringsAsFactors = FALSE
)
# 4. Convert the matrix to long format, join with our hierarchy map descriptions
functional_long <- as.data.frame(t(clean_numeric_matrix)) %>%
  tibble::rownames_to_column("KO_Number") %>%
  tidyr::pivot_longer(cols = -KO_Number, names_to = "SampleID", values_to = "Count") %>%
  dplyr::left_join(hierarchy_map, by = "KO_Number")

# 5. Create a descriptive, unique column name combining Level information and KO ID
# This prevents R from creating duplicate column headers during the widening step
functional_long <- functional_long %>%
  dplyr::mutate(Pathway_Column_Header = paste(Level1, KO_Number, sep = "_"))

# 6. Pivot the table wide so that Rows = Samples, and Columns = Your specific functions!
functional_wide_metadata <- functional_long %>%
  dplyr::select(SampleID, Pathway_Column_Header, Count) %>%
  tidyr::pivot_wider(names_from = Pathway_Column_Header, values_from = Count) %>%
  as.data.frame()
rownames(functional_wide_metadata) <- functional_wide_metadata$SampleID


# now we have the functional table in a interpretable data format, for use in further analysis. now we aim to combine it with our 
# metadata, so we can attach it in the end and use it in our visualisation and data analysis.

# 1. Pull the original metadata dataframe from your phyloseq object
metadata_df          <- as(sample_data(ps_work), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# 2. Extract only the first 10 columns for easy viewing
# We use unique() to make sure 'SampleID' is included even if it wasn't in the first 10 columns originally
columns_to_keep         <- unique(c(colnames(metadata_df)[1:100], "SampleID"))
metadata_pre_copy_paste <- metadata_df[, columns_to_keep, drop = FALSE]
#Fix the sample order before merging:

# 3. Merge the compact metadata with your functional columns, in the right column using the left join function
metadata_to_copy_paste <- metadata_pre_copy_paste %>%
  dplyr::left_join(functional_wide_metadata, by = "SampleID")

# Optional: Put SampleID as the very first column for easy Excel copy-pasting
metadata_to_copy_paste <- metadata_to_copy_paste[, c("SampleID", setdiff(colnames(metadata_to_copy_paste), "SampleID"))]

#now we fix the lexicographical sorting which R uses as a baseline if it does not see the labels as numeric:
# Convert the order column to numeric and arrange the data frame
metadata_sorted <- metadata_to_copy_paste[order(as.numeric(metadata_to_copy_paste$order_treatment)), ]

#Before saving the object, fix the collumn names, by replacing all spaces and special characters with an _ 
colnames(metadata_sorted) <- gsub("/", "_", colnames(metadata_sorted))
colnames(metadata_sorted) <- gsub(",", "", colnames(metadata_sorted))
colnames(metadata_sorted) <- gsub(" ", "_", colnames(metadata_sorted))

head(metadata_sorted)
# 4. Create an output directory if it doesn't exist, and write out the file
if(!dir.exists("results/all_wetlands/")) {
  dir.create("results/all_wetlands/", recursive = TRUE)
}

# 5. Write the tab-separated file safely (Using the correct backslash \t)
write.table(metadata_sorted, 
            sep = "\t", dec = ",",
            file = "results/all_wetlands/Extended_Functional_Respiration_Metadata.tsv", 
            row.names = FALSE,
            quote = TRUE)

message("Resulting table which can be copy pasted can be found as Extended_Functional_Respiration_Metadata.tsv in the results/all_wetlands/...")
