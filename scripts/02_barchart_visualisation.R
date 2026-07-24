# 08_barchart_visualisation.R

#' Prepare data for Relative Abundance Plots (Single or Dual Level Mode)
prepare_abundance_data <- function(ps, high_level = "phylum", low_level = "order", 
                                   facet_var = "Sample", threshold = 1) {
  library(dplyr)
  
  # 1. Transform and Melt
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
  df <- psmelt(ps_rel)
  
  # Dynamic Check FIX: Clean evaluation that yields a single TRUE or FALSE
  has_low_level <- !is.null(low_level) && !any(is.na(low_level)) && low_level != ""
  
  # 2. Average by group if specified
  if (facet_var != "Sample") {
    group_vars <- c("OTU", facet_var, high_level)
    if (has_low_level) {
      group_vars <- c(group_vars, low_level)
    }
    
    df <- df %>%
      dplyr::group_by(across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(Abundance = mean(Abundance), .groups = "drop")
  }
  
  # 3. CRITICAL FIXED LOGIC: Calculate full taxonomic group totals BEFORE threshold check
  # Otherwise, individual small OTUs will fail the threshold and everything becomes "Others"
  group_total_vars <- c(facet_var, high_level)
  if (has_low_level) {
    group_total_vars <- c(group_total_vars, low_level)
  }
  
  df <- df %>%
    dplyr::group_by(across(dplyr::all_of(group_total_vars))) %>%
    dplyr::mutate(Group_Total_Abund = sum(Abundance)) %>%
    dplyr::ungroup()
  
  # 4. Define Labels using the true Group totals
  df <- df %>%
    dplyr::mutate(
      Level_High = as.character(!!sym(high_level)),
      Plot_High = ifelse(Group_Total_Abund >= threshold, Level_High, "Others")
    )
  
  # Inject conditional label assignment path
  if (has_low_level) {
    df <- df %>%
      dplyr::mutate(
        Level_Low  = as.character(!!sym(low_level)),
        Plot_Low   = ifelse(Group_Total_Abund >= threshold, Level_Low, "Minor Taxa"),
        Hierarchical_Label = ifelse(Plot_High == "Others", "Others", paste(Plot_High, Plot_Low, sep = " - "))
      )
  } else {
    df <- df %>%
      dplyr::mutate(
        Plot_Low   = "Single_Level",
        Hierarchical_Label = ifelse(Plot_High == "Others", "Others", Plot_High)
      )
  }
  
  # Section to put the others on the bottom of the graph
  df <- df %>%
    dplyr::group_by(!!sym(facet_var), Plot_High, Hierarchical_Label) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
  
  all_labels <- unique(df$Hierarchical_Label)
  major_taxa <- sort(setdiff(all_labels, "Others"))
  target_levels <- c(major_taxa, "Others")
  
  df$Hierarchical_Label <- factor(df$Hierarchical_Label, levels = target_levels)
  
  return(df)
}

#' Create shades for a base color (Kept exactly as yours)
create_shades <- function(base_color, n) {
  if (n <= 1) return(base_color)
  grad_pal <- scales::seq_gradient_pal("#F0F0F0", base_color, "Lab")
  return(grad_pal(seq(0.3, 1, length.out = n)))
}

#' Generate Shaded Color Palette (Kept exactly as yours)
get_shaded_palette <- function(df, pal_name = "ggthemes::Tableau_10") {
  library(paletteer)
  
  high_groups <- df %>% 
    dplyr::filter(Plot_High != "Others") %>% 
    dplyr::pull(Plot_High) %>% 
    unique() %>% 
    sort()
  
  n_groups <- length(high_groups)
  if (n_groups == 0) {
    stop("Error: No taxa passed the abundance threshold filter. Try reducing your threshold in the config.")
  }
  
  raw_pal <- paletteer_d(pal_name)
  interp_pal <- grDevices::colorRampPalette(raw_pal)(n_groups)
  names(interp_pal) <- high_groups
  
  unique_labels <- df %>%
    dplyr::select(Plot_High, Hierarchical_Label) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Plot_High == "Others", Plot_High, Hierarchical_Label)
  
  color_map <- unique_labels %>%
    dplyr::group_by(Plot_High) %>%
    dplyr::mutate(
      Shade = if (unique(Plot_High) == "Others") {
        "#A9A9A9" 
      } else {
        create_shades(interp_pal[unique(Plot_High)], dplyr::n())
      }
    ) %>%
    dplyr::ungroup()
  
  return(setNames(color_map$Shade, color_map$Hierarchical_Label))
}