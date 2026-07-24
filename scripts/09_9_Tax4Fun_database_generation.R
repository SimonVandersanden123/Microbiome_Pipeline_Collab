#Script to pull WGS from the mibirem whole genome database for generation of a reference database for functional assignment of the taxa found in the sites
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# This is something that could be done in the future to have a mibirem specific database for functional assignment of the microbiomes 
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
library(biomartr)
library(jsonlite)
library(R.utils)
library(Tax4Fun2)


###########
####I first downloaded the reference genomes from the NCBI dataset using this command in an anaconda prompt 
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA629478
#something like this :# 1. Move into the folder where you downloaded the file
#  cd C:\ncbi_work
#
## 2. Run the command to get the dehydrated skeleton package
#.\datasets download genome accession PRJNA629478 --dehydrated --filename midas_dehydrated.zip
#
## 3. Extract the skeleton structure
#Expand-Archive -Path .\midas_dehydrated.zip -DestinationPath .\midas_package
#
## 4. Trigger the actual download of the 1,000+ genomes
#.\datasets rehydrate --directory .\midas_package\
###########
# Now we try to use these downloaded genomes to train a new reference database
# the reference genomes are in here:work from here
#C:\ncbi_work\midas_package\ncbi_dataset\data

# Define target input and output files
assembly_file <- "PRJEB88388_AssemblyDetails.txt"
output_dir    <- "mibirem_master_genomes"

# Create the folder where genomes will be stored
dir.create(output_dir, showWarnings = FALSE)
################################################################################
# 2. FILTER & EXTRACT ACCESSIONS (Replaces 'awk')
################################################################################
cat("Reading and filtering assembly details...\n")

# Read the raw lines from the text file
lines <- readLines(assembly_file)

# Find rows that start with a GCA accession number
gca_lines <- lines[grep("^GCA_", lines)]

# Clean up tabs/spaces to separate the columns accurately
parsed_data <- do.call(rbind, strsplit(gca_lines, "\\t+"))

# Extract columns (Column 1 = Accession, Column 2 = Assembly Level)
accessions <- parsed_data[, 1]
levels     <- parsed_data[, 2]

# OPTION A: If you want ALL genomes (Complete + Contigs):
target_accessions <- unique(accessions)

# OPTION B: If you want ONLY Complete Genomes (Recommended by Claude for cleaner DB):
# target_accessions <- unique(accessions[levels == "Complete genome"])

cat("Found", length(target_accessions), "target accessions to download.\n")


#Check if all the desired genomes are there

list.files("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab/tax4fun2/Tax4Fun2_Custom_Reference_Database/mibirem_master_genomes")
################################################################################
# Now all the genomes are downloaded and put in the map 
################################################################################
# https://github.com/songweizhi/Tax4Fun2_short_tutorial/blob/master/db_default_user.md 
# 1. Set Working Directory
setwd("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab")
pwd_user_data  <- "tax4fun2/Tax4Fun2_Custom_Reference_Database/mibirem_master_genomes" # The folder with your .fna files
pwd_ref_data   <- "tax4fun2/Tax4Fun2_ReferenceData_v2"   # path to the reference data
name_user_data <- "tax4fun2/Tax4Fun2_Custom_Reference_Database/mibirem_ref_db"         # What you want to name your custom DB


# 1. Normalize your paths so Windows command-line binaries can read them
pwd_user_data_clean <- normalizePath("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab/tax4fun2/Tax4Fun2_Custom_Reference_Database/mibirem_master_genomes")
pwd_ref_data_clean  <- normalizePath("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab/tax4fun2/Tax4Fun2_ReferenceData_v2")


################################################################################
# STEP A: Extract 16S (SSU) rRNA sequences from your downloaded genomes
################################################################################
cat("Extracting 16S rRNA sequences from genomes...\n")

extractSSU(
  genome_folder = pwd_user_data, 
  file_extension = "fna", 
  path_to_reference_data = pwd_ref_data
)
From this point on the reference documentation does not work
It is just suppose to be this short
library(Tax4Fun2)

setwd("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab")
pwd_user_data  <- "tax4fun2/Tax4Fun2_Custom_Reference_Database/mibirem_master_genomes" # The folder with your .fna files
pwd_ref_data   <- "tax4fun2/Tax4Fun2_ReferenceData_v2"   # path to the reference data
name_user_data <- "tax4fun2/Tax4Fun2_Custom_Reference_Database/mibirem_ref_db"         # What you want to name your custom DB

assignFunction(genome_folder = pwd_user_data, file_extension = "fna", path_to_reference_data = pwd_ref_data, num_of_threads = 1, fast = TRUE)


# modify the following 4 lines as needed
pwd_ref_data   = 'Tax4Fun2_ReferenceData_v2'   # path to Tax4Fun2's default database, need to be decompressed before use
pwd_user_data  = 'genome_folder'               # path to the folder that holds the reference genome/MAG files
name_user_data = 'name_of_user_database'       # specify the name of the generated database, specify only the name, do not include path here!
gnm_ext        = 'fna'                         # extension of the genome files

# Generate your own database
extractSSU(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data)

generateUserData(path_to_reference_data = pwd_ref_data, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data, SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")

