
_______________TAx4FUN2____________________
Installation and setup
```{r}
#Run these steps for the installation of the package and downloading the reference databases:
install.packages(pkgs = "C:/Users/lucp12540/Downloads/Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)
# 1. Point R to your exact pipeline folder (using forward slashes)
setwd("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab/tax4fun2")
# 2. Extract the reference data archive directly into this folder
untar("C:/Users/lucp12540/Downloads/Tax4Fun2_ReferenceData_v2.tar.gz", exdir = ".")

# 3. Load the package
library(Tax4Fun2)
# 1. Point R to your exact pipeline folder (using forward slashes)
setwd("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab/tax4fun2")
# Force the package to download the official native Windows BLAST binaries
buildDependencies(
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
  install_suggested_packages = TRUE, 
  use_force = TRUE
)

# 1. Point R to your exact active workspace folder
setwd("C:/Users/lucp12540/Documents/GitHub/Microbiome_Pipeline_Collab")

# 3. Rename the extracted folder to 'blast_bin' so Tax4Fun2 can find 'blastn.exe'
file.rename(
  from = "tax4fun2/Tax4Fun2_ReferenceData_v2/ncbi-blast-2.9.0+", 
  to = "tax4fun2/Tax4Fun2_ReferenceData_v2/blast_bin"
)

```