# Tax4Fun object generation "09_Tax4Fun_object_generation.R"
#load dependencies
library(Tax4Fun2)
# 1. Extract the unique sequences (the column names of your OTU table)
asv_sequences <- colnames(otu_table(ps_work))
# 2. Write the Fasta file (Sequence ID and Sequence are identical here)
fasta_lines <- list()
for (i in 1:length(asv_sequences)) {
  fasta_lines[[paste0(">Seq_", i)]] <- paste0(">Seq_", i, "\n", asv_sequences[i])
}
writeLines(unlist(fasta_lines), "tax4fun2/mibi_pipeline_otus.fasta")
# 3. Extract the Abundance Table
# Tax4Fun2 requires Sample IDs as rows/columns matching the orientation. 
# Since your data has taxa_are_rows = FALSE, we need to transpose it so sequences are rows.
otu_matrix <- t(as.matrix(otu_table(ps_work)))
# Replace the actual long sequences with the sequence keys ("Seq_1", "Seq_2"...) to match the Fasta file
rownames(otu_matrix) <- paste0("Seq_", 1:length(asv_sequences))
# Save out as a tab-separated file
write.table(otu_matrix, "tax4fun2/mibi_pipeline_otu_table.txt", sep="\t", col.names=NA, quote=FALSE)
# Now the objects are generated to be put into the specified functions:

# predict functions
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)