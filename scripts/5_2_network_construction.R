# =====================================================================
# 5.2 Microbial Network Construction Script
# =====================================================================
#
# Network construction step, there are multiple ways of constructing the network: 
# =====================================================================
# There are 3 main ways of associating the taxa of interest:
# Correlation, proportionality and conditional dependence
# From each of these measures, there are also different ways of data filtering, zero handling and normalization methods.
# I have included one section for each one of these major association methods to generate the networks, since the choice of association has a great impact on the structure of the final network.
#
# Specific parameters customisation is mentioned in the 'help' files associated with the function, for some important parameters, short information is included in the current script.
#                                         # help("netConstruct")
# An important one which is important for the visualisation: dissFunc: there are multple measures to calculate the dissimilarity between the different taxa which will determine the position of the taxa within the final visualisation. And therefore also the interpretation or intuition of the resulting network. 
#   #-unsigned: leading to low distance between strongly associated taxa (both positively associated as well as negatively associated) 
#   #-signed: where the distance between the nodes is highest for strongly negative associated taxa. (default)
#   #-SignedPos, signed distance with setting negative assocations to zero.
#   #-There is also TOMdiss: dissimilarity measure available based on the topological overlap matrix (TOM), instead of looking at the relationships between the taxa in isolation, it evaluates the shared neighborhood context. (if taxa share a lot of the same neighbors, they will be closer to one another)
#
# =====================================================================
# Theoretical Background for Association Methods
# =====================================================================
# 1. Correlation based methods (SparCC):
# We start off with Correlation, since this is the most intuitive and 'straight forward method:
#   ## Correlation-based networks represent co-occurrence patterns (marginal associations) rather than direct biological interactions. In these networks, if Taxon A and Taxon B both thrive in the same environmental niche (e.g., sharing a specific depth or pH preference), they will show a strong positive correlation. However, this correlation is often a "ghost" or spurious association driven by a shared environmental confounder rather than a direct relationship.
#
#   ## To identify actual direct relationships, conditional dependence methods (such as SPRING) are used. These methods evaluate whether Taxon A and Taxon B remain associated after conditioning on (controlling for) all other taxa in the dataset (Taxon C,D,...). This statistical control actively filters out indirect paths and shared environmental effects.
#
# ## Network Construction using SparCC: For our correlation-based analysis, we utilized SparCC (Sparse Correlations for Compositional data). Traditional correlation metrics (like Pearson or Spearman) yield spurious results on microbiome datasets due to their compositional nature (where relative abundances must sum to 1). SparCC is specifically designed to be compositionally aware, utilizing log-ratio transformations to estimate true correlation coefficients without compositionality bias.
# ### Here normalisation and zero replacement methods are automatically selected so there is no need to specify it 
#
# ---------------------------------------------------------------------
# 2. Conditional dependence based method (SPRING):
# ## Conditional Dependence (Direct Relationships): This approach determines if Taxon A and Taxon B are still connected once we statistically account for the rest of the community (Taxon C). It strips away these indirect, environmentally-driven "ghost" connections to reveal the true backbone of the network.
#
# ## We do the network constuction using the SPRING, Unlike traditional correlation metrics that yield spurious results on compositional microbiome data (where relative abundances must sum to 1), SPRING is fully compositionally aware.This allows SPRING to estimate a sparse precision matrix representing true conditional independence, successfully distinguishing direct ecological interactions from indirect background noise.
# ### Here normalisation and zero replacement methods are automatically selected so there is no need to specify it 
#
# ---------------------------------------------------------------------
# 3. Proportionality based method (propr / rho):
# ## Network Construction using Proportionality (rho)
# We construct the microbial networks using **Proportionality (rho)**. Unlike traditional correlation metrics (such as Pearson or Spearman) that yield spurious (fake) associations on compositional relative abundance data, the proportionality measure rho is fully compositionally aware and sub-compositionally coherent. It is a symmetric measure with values scaled between -1 and +1 (where +1 represents perfect positive proportionality), making it directly comparable to standard correlation coefficients and suitable for standard network topology algorithms.
#
# #### Connection to Pearson Correlation
# The formula for rho is highly analogous to the classic **Pearson Correlation Coefficient (r)**. 
# Recall that Pearson's correlation scales covariance using the geometric mean of the variances:
# r(Ai, Aj) = cov(Ai, Aj) / square_root([var(Ai)]^2 * [var(Aj)]^2)
# In contrast, proportionality rho scales the relationship using the arithmetic mean of the variances. It is calculated as: 
# rho(log(x_i),log(x_j) = 1 - [ var( log(x_i / x_j) ) / ( var(log(x_i)) + var(log(x_j)) ) ]
#
# #### Breaking Down the Formula Components:
# * **log(x_i / x_j)**: Takes the logarithm of the ratio of the two taxa. This is mathematically equivalent to the difference of their logs: log(x_i) - log(x_j). if this log(x_i)- log(y_i) = a constant value, this means that the ratio's of these taxa are roughly equal overall and, then the variance of this term is low.=>
# * **var( log(x_i / x_j) )**: Calculates the Log-Ratio Variance (VLR). If the ratio between two taxa stays highly consistent across different samples/conditions, this variance will be very close to 0.
# * **var(log(x_i)) + var(log(x_j))**: Represents the sum of the individual variances of the log-transformed data, acting as a scaling/standardising factor. (Note: In practice, to make this compositionally safe, NetCoMi calculates this using Centered Log-Ratio (CLR) transformed values).
#
# #### The Three Core Scenarios for interpreting the rho scores:
# ##### 1. Perfect Positive Proportionality (rho = 1)
# If Taxon x is always for example exactly 3 times more abundant than Taxon y, their ratio x/y is a constant (3). 
# * The log of a constant is a constant (log(3) is approximately 1.099).
# * The variance of a constant across samples is 0. Plugging this into the numerator of the formula yields:
# rho = 1 - [ 0 / (var(log(x)) + var(log(y))) ] = 1 - 0 = 1
#
# ##### 2. Complete Independence (rho = 0)
# If the two genera have absolutely nothing to do with each other, the variance of their difference is simply the sum of their individual variances: 
# var(log(x) - log(y)) = var(log(x)) + var(log(y))
# * Under this condition, the numerator becomes identical to the denominator.Therefore the fraction simplifies to 1, yielding:
# rho = 1 - 1 = 0
#
# ##### 3. Perfect Inverse Proportionality (rho = -1)
# If one taxon increases while the other decreases, their covariance is highly negative. The variance of their difference expands because of the algebraic rule:
# var(log(x) - log(y)) = var(log(x)) + var(log(y)) - [ 2 * cov(log(x), log(y)) ]
# * Because covariance is negative, subtracting it becomes addition: - [ 2 * (-covariance) ] = + [ 2 * absolute_value_of_covariance ]. 
# * This makes the variance of their difference (the numerator) significantly larger than the sum of their individual variances (the denominator).
# * Due to mathematical boundaries, the numerator can never exceed twice the value of the denominator. Thus, the fraction is capped at 2:
# rho = 1 - 2 = -1
# * This lower limit of -1 represents a perfect reciprocal relationship. Because the raw data is first CLR-transformed, even exponential differences between competing bacteria are successfully projected into this symmetric, linear -1 to +1 space.
#
# ### To achieve this, rho modifies the log-ratio variance using Centered Log-Ratio (clr) transformed data. This provides a scaled, robust association metric that is mathematically immune to compositionality bias, allowing us to reconstruct true co-abundance patterns without the risk of false-positive correlations.
# =====================================================================

library(NetCoMi)
library(phyloseq)
library(remotes)
library(propr)

# ---------------------------------------------------------------------
# Step 1. Prepare Grouping Factor
# ---------------------------------------------------------------------
ps_net_input <- net_build_params$ps_input

group_vector <- NULL
if (!is.null(net_build_params$group_var) && net_build_params$group_var %in% colnames(sample_data(ps_net_input))) {
  group_vector <- as.character(sample_data(ps_net_input)[[net_build_params$group_var]])
  message("Grouping variable set to: ", net_build_params$group_var)
} else {
  message("No grouping variable specified or found. Constructing single network.")
}

# ---------------------------------------------------------------------
# Step 2. Execute Selected Association Method
# ---------------------------------------------------------------------
selected_method <- tolower(net_build_params$method)
message("Building network using association method: ", selected_method)

if (selected_method %in% c("sparcc", "correlation")) {
  
  # --- SparCC (Correlation) ---
  network_to_analyse <- netConstruct(
    data        = ps_net_input,
    measure     = "sparcc",
    taxRank     = net_build_params$tax_level,
    filtTax     = net_build_params$filter_tax_method,
    filtTaxPar  = list(highestFreq = net_build_params$top_taxa_n),
    group       = group_vector,
    sparsMethod = net_build_params$spars_method,
    thresh      = net_build_params$thresh_sparcc,
    dissFunc    = net_build_params$diss_func,
    verbose     = net_build_params$verbose,
    seed        = net_build_params$seed
  )
  
} else if (selected_method %in% c("spring", "conditional", "conditional_dependence")) {
  
  # --- SPRING (Conditional Dependence) ---
  network_to_analyse <- netConstruct(
    data        = ps_net_input,
    measure     = "spring",
    taxRank     = net_build_params$tax_level,
    filtTax     = net_build_params$filter_tax_method,
    filtTaxPar  = list(highestFreq = net_build_params$top_taxa_n),
    filtSamp    = "none",
    group       = group_vector,
    measurePar  = list(
      nlambda = net_build_params$spring_nlambda,
      rep.num = net_build_params$spring_rep_num,
      Rmethod = net_build_params$spring_rmethod
    ),
    dissFunc    = net_build_params$diss_func,
    verbose     = net_build_params$verbose,
    seed        = net_build_params$seed
  )
  
} else if (selected_method %in% c("propr", "proportionality")) {
  
  # --- Proportionality (rho) ---
  network_to_analyse <- netConstruct(
    data        = ps_net_input,
    measure     = "propr",
    taxRank     = net_build_params$tax_level,
    filtTax     = net_build_params$filter_tax_method,
    filtTaxPar  = list(highestFreq = net_build_params$top_taxa_n),
    filtSamp    = "none",
    group       = group_vector,
    sparsMethod = net_build_params$spars_method,
    thresh      = net_build_params$thresh_propr,
    dissFunc    = net_build_params$diss_func,
    verbose     = net_build_params$verbose,
    seed        = net_build_params$seed
  )
  
} else {
  stop("Invalid method specified in net_build_params$method! Choose 'sparcc', 'spring', or 'propr'.")
}

# ---------------------------------------------------------------------
# Step 3. Standardized Output Assignment & Caching
# ---------------------------------------------------------------------

# Expose output object into workspace
assign("network_to_analyse", network_to_analyse, envir = .GlobalEnv)

if (isTRUE(net_build_params$save_output) && !is.null(net_build_params$output_rds_path)) {
  dir.create(dirname(net_build_params$output_rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(network_to_analyse, file = net_build_params$output_rds_path)
  message("Saved constructed network object to: ", net_build_params$output_rds_path)
}

message("Network construction complete. Output object ready: 'network_to_analyse'")