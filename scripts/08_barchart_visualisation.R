# 08_barchart_visualisation.R
#' Prepare data for Relative Abundance Plots
prepare_abundance_data <- function(ps, high_level = "phylum", low_level = "order", 
                                   facet_var = "Sample", threshold = 1) {
  
  # 1. Transform and Melt
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
  df <- psmelt(ps_rel)
  
  # 2. Average by group if specified
  if (facet_var != "Sample") {
    df <- df %>%
      dplyr::group_by(OTU, !!sym(facet_var), !!sym(high_level), !!sym(low_level)) %>%
      dplyr::summarise(Abundance = mean(Abundance), .groups = "drop")
  }
  # 3. Define Labels and put taxa which do not meet the threshold into 'Others'
  df <- df %>%
    dplyr::mutate(
      Level_High = as.character(!!sym(high_level)),
      Level_Low  = as.character(!!sym(low_level)),
      Plot_High = ifelse(Abundance >= threshold, Level_High, "Others"),
      Plot_Low  = ifelse(Abundance >= threshold, Level_Low, "Minor Taxa"),
      Hierarchical_Label = ifelse(Plot_High == "Others", "Others", paste(Plot_High, Plot_Low, sep = " - "))
    )
  #Section to put the others on the bottom of the graph:
  df <- df %>%
    dplyr::group_by(!!sym(facet_var), Plot_High, Hierarchical_Label) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
  all_labels <- unique(df$Hierarchical_Label)
  # Identify all labels that ARE NOT "Others" and sort them alphabetically
  major_taxa <- sort(setdiff(all_labels, "Others"))
  target_levels <- c(major_taxa ,"Others")
  # Combine: "Others" first (bottom of stack), followed by major taxa
  df$Hierarchical_Label <- factor(df$Hierarchical_Label, levels = target_levels)
  # ---------------------------
  return(df)
}

#' Create shades for a base color
create_shades <- function(base_color, n) {
  if (n <= 1) return(base_color)
  # Gradient from nearly white (0.3) to the full base color (1)
  grad_pal <- scales::seq_gradient_pal("#F0F0F0", base_color, "Lab")
  return(grad_pal(seq(0.3, 1, length.out = n)))
}

#' Generate Shaded Color Palette.
get_shaded_palette <- function(df, pal_name = "ggthemes::Tableau_10") {
  
  # 1. Get unique High-Level groups
  high_groups <- df %>% 
    dplyr::filter(Plot_High != "Others") %>% 
    dplyr::pull(Plot_High) %>% 
    unique() %>% 
    sort()
  
  n_groups <- length(high_groups)
  # 2. Handle Palette Overflow via Interpolation
  # We extract the palette and use colorRampPalette to stretch it to n_groups
  raw_pal <- paletteer_d(pal_name)
  interp_pal <- grDevices::colorRampPalette(raw_pal)(n_groups)
  names(interp_pal) <- high_groups
  
  # 3. Create the Shade Mapping
  # Use dplyr::select explicitly to avoid namespace errors
  unique_labels <- df %>%
    dplyr::select(Plot_High, Hierarchical_Label) %>%
    dplyr::distinct() %>%
    # Sort so 'Others' is handled separately or last
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