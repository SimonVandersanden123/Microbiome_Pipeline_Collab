# 03_alpha_plots.R
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Preparation
plot_data <- alpha_output %>%
  mutate(across(all_of(c(plot_x_axis, plot_color_by)), as.factor)) %>%
  pivot_longer(cols = all_of(plot_metrics), 
               names_to = "Metric", 
               values_to = "Value")

# 2. Dynamic Plotting
alpha_boxplot <- ggplot(plot_data, aes(x = .data[[plot_x_axis]], y = Value, fill = .data[[plot_color_by]])) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), 
             size = 1.5, alpha = 0.5, color = "black") +
  facet_wrap(~Metric, scales = "free_y", ncol = plot_cols) +
  theme_bw() +
  labs(title = "Microbial Community Trait Analysis",
       x = plot_x_axis, y = "Index Value", fill = plot_color_by) +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

message("Alpha plot generated.")