# filter for edge values only for CM only
EE_plate_edges <- EE_plate %>% 
  filter(row %in% c("A", "B", "P", "O")) %>% 
  filter(media_type %in% "CM")

# A VS P
# filter for A
EE_plate_edge_A <- EE_plate_edges %>% 
  filter(row %in% "A")

# filter for P
EE_plate_edge_P <- EE_plate_edges %>% 
  filter(row %in% "P")

# A: filter for compound + area + 7.5K
EE_7.5K_area_A <- EE_plate_edge_A %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) 
  # select(media_type, lps, compound, area7_5k_nuc)

EE_7.5K_area_A_plot <- EE_7.5K_area_A %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  labs(title = "Min/Max Data (Area) from 'A' Wells using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Key: Treatment") +
  theme_minimal() +
  theme(legend.position = "none")


# P: filter for compound + area + 7.5K
EE_7.5K_area_P <- EE_plate_edge_P%>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) 
# select(media_type, lps, compound, area7_5k_nuc)

EE_7.5K_area_P_plot <- EE_7.5K_area_P %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  labs(title = "Min/Max Data (Area) from 'P' Wells using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Key: Treatment") +
  theme_minimal() +
  theme(legend.position = "none")


# Find the maximum value across both plots
max_value <- max(c(max(EE_7.5K_area_A$area7_5k_nuc), max(EE_7.5K_area_P$area7_5k_nuc)))

# Set the same y-axis limits for both plots
EE_7.5K_area_A_plot <- EE_7.5K_area_A_plot + coord_cartesian(ylim = c(0, max_value))
EE_7.5K_area_P_plot <- EE_7.5K_area_P_plot + coord_cartesian(ylim = c(0, max_value))

library(gridExtra)
grid.arrange(EE_7.5K_area_A_plot, EE_7.5K_area_P_plot, ncol=2)




# B VS O
# filter for B
EE_plate_edge_B <- EE_plate_edges %>% 
  filter(row %in% "B")

# filter for O
EE_plate_edge_O <- EE_plate_edges %>% 
  filter(row %in% "O")

# B: filter for compound + area + 7.5K
EE_7.5K_area_B <- EE_plate_edge_B %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) 
# select(media_type, lps, compound, area7_5k_nuc)

EE_7.5K_area_B_plot <- EE_7.5K_area_B %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  labs(title = "Min/Max Data (Area) from 'B' Wells using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Key: Treatment") +
  theme_minimal() +
  theme(legend.position = "none")


# O: filter for compound + area + 7.5K
EE_7.5K_area_O <- EE_plate_edge_O %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) 
# select(media_type, lps, compound, area7_5k_nuc)

EE_7.5K_area_O_plot <- EE_7.5K_area_O %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  labs(title = "Min/Max Data (Area) from 'O' Wells using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Key: Treatment") +
  theme_minimal() +
  theme(legend.position = "none")


# Find the maximum value across both plots
max_value <- max(c(max(EE_7.5K_area_B$area7_5k_nuc), max(EE_7.5K_area_O$area7_5k_nuc)))

# Set the same y-axis limits for both plots
EE_7.5K_area_B_plot <- EE_7.5K_area_B_plot + coord_cartesian(ylim = c(0, max_value))
EE_7.5K_area_O_plot <- EE_7.5K_area_O_plot + coord_cartesian(ylim = c(0, max_value))

library(gridExtra)
grid.arrange(EE_7.5K_area_B_plot, EE_7.5K_area_O_plot, ncol=2)
