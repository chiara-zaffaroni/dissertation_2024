# filter for CM (edgeless)
# filter for area and 7.5K 
# max signal (plateau): LPS/dmso
# min TAK signal: LPS/TAK3
# min noLPS signal: noLPS/dmso

# MIN/MAX EDGELESS 7.5K
# filter for CM + area + 7.5K
EE_7.5K_non_edges <- EE_plate_edgefree %>% 
  filter(media_type %in% "CM") %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) %>% 
  # remove outliers
  filter(area7_5k_nuc < 300) %>% 
  filter(area7_5k_nuc > 10)
  # select(media_type, lps, compound, area7_5k_nuc)

# non-edge plot 
EE_7.5K_non_edges_plot <- EE_7.5K_non_edges %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "Min/Max Data (Area) from Non-Edge Wells using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Treatment") +
  theme_minimal() +
  theme(legend.position = "none")

# Group by combinations and calculate summary statistics
summary_EE_7.5K_area <- EE_7.5K_non_edges %>%
  group_by(lps, compound, .groups = "drop") %>%
  summarise(
    N = n(),
    Mean = mean(area7_5k_nuc),
    Median = median(area7_5k_nuc),
    Min = min(area7_5k_nuc),
    Max = max(area7_5k_nuc),
    SD = sd(area7_5k_nuc),
    SEM = sd(area7_5k_nuc) / sqrt(length(area7_5k_nuc)),
    CV = sd(area7_5k_nuc) / mean(area7_5k_nuc) * 100
  ) %>% 
  gt()

# outliers
install.packages("robustbase")
library(robustbase)
model <- lmrob(area7_5k_nuc ~ interaction(lps, compound), data = EE_7.5K_area)
your_data_filtered <- EE_7.5K_area[model$rweights > 0, ]

# MIN/MAX WITH EDGES 7.5K
# filter for edge values only for CM only for 7.5K
EE_7.5K_edges <- EE_plate %>% 
  filter(row %in% c("A", "B", "P", "O")) %>% 
  filter(media_type %in% "CM") %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) %>% 
  # remove outlier
  filter(area7_5k_nuc < 120)  
 
# edge plot
EE_7.5K_edges_plot <- EE_7.5K_edges %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "Min/Max Data (Area) from Edge Wells using 7.5K threshold",
       x = "Treatment",
       fill = "Treatment") +
  theme_minimal() +
  theme(axis.title.y=element_blank())


# Find the maximum value across both plots
max_value <- max(c(max(EE_7.5K_non_edges$area7_5k_nuc), max(EE_7.5K_edges$area7_5k_nuc)))

# Set the same y-axis limits for both plots
EE_7.5K_non_edges_plot <- EE_7.5K_non_edges_plot + coord_cartesian(ylim = c(0, max_value))
EE_7.5K_edges_plot <- EE_7.5K_edges_plot + coord_cartesian(ylim = c(0, max_value))

grid.arrange(EE_7.5K_non_edges_plot, EE_7.5K_edges_plot, ncol=2)

# FOR EACH DRUG