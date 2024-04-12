# filter for edge values only for CM only
EE_plate_edges <- EE_plate %>% 
  filter(row %in% c("A", "B", "P", "O")) %>% 
  filter(media_type %in% "CM")

# filter for A+P
EE_plate_edges_AP <- EE_plate_edges %>% 
  filter(row %in% c("A", "P"))

# filter for B+O
EE_plate_edges_BO <- EE_plate_edges %>% 
  filter(row %in% c("B", "O"))

# AP: filter for compound + area + 7.5K
EE_7.5K_area_AP <- EE_plate_edges_AP %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) %>% 
  select(media_type, lps, compound, area7_5k_nuc)

EE_7.5K_area_AP_plot <- EE_7.5K_area_AP %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  labs(title = "Min/Max Data (Area) from Outermost Edge Wells (A+P) using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Key: Treatment") +
  theme_minimal() +
  theme(legend.position = "none")

# BO: filter for compound + area + 7.5K
EE_7.5K_area_BO <- EE_plate_edges_BO %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(area7_5k_nuc = as.numeric(area7_5k_nuc)) %>% 
  filter(area7_5k_nuc < 120) %>% 
  select(media_type, lps, compound, area7_5k_nuc)


EE_7.5K_area_BO_plot <- EE_7.5K_area_BO %>% 
  ggplot(aes(x = interaction(lps, compound), y = area7_5k_nuc, fill = interaction(lps, compound))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  labs(title = "Min/Max Data (Area) from Second Outermost Edge Wells (B+O) using 7.5K threshold",
       x = "Treatment",
       y = "CD38 area/nuclei",
       fill = "Key: Treatment") +
  theme_minimal()

# Find the maximum value across both plots
max_value <- max(c(max(EE_7.5K_area_AP$area7_5k_nuc), max(EE_7.5K_area_BO$area7_5k_nuc)))

# Set the same y-axis limits for both plots
EE_7.5K_area_AP_plot <- EE_7.5K_area_AP_plot + coord_cartesian(ylim = c(0, max_value))
EE_7.5K_area_BO_plot <- EE_7.5K_area_BO_plot + coord_cartesian(ylim = c(0, max_value))

library(gridExtra)
grid.arrange(EE_7.5K_area_AP_plot, EE_7.5K_area_BO_plot, ncol=2)



# AP SUMMARY STATS
summary_EE_7.5K_area_AP <- EE_7.5K_area_AP %>%
  group_by(lps, compound, Treatment = " ") %>%
  summarise(
    N = n(),
    Mean = mean(area7_5k_nuc),
    Median = median(area7_5k_nuc),
    # Min = min(area7_5k_nuc),
    # Max = max(area7_5k_nuc),
    SD = sd(area7_5k_nuc),
    SEM = sd(area7_5k_nuc) / sqrt(length(area7_5k_nuc)),
    CV = sd(area7_5k_nuc) / mean(area7_5k_nuc) * 100
  ) %>% 
  gt(rowname_col = "Treatment") %>% 
  tab_header(title = "Analysis: Outermost (A+P) Edges") %>%
  fmt_number(decimals = 2) %>% 
  tab_style(
    style = cell_text(
      size = "smaller",
      weight = "bold",
    ),
    locations = cells_body(columns = "CV")
  ) %>% 
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = "CV")
  )

# 
#   tab_style(style = list(
#     cell_fill(color = "lightcyan"),
#     cell_text(weight = "bold")
#   ),
#   locations = cells_body(
#     columns = CV
#   ))

# BO SUMMARY STATS
summary_EE_7.5K_area_BO <- EE_7.5K_area_BO %>%
  group_by(lps, compound, Treatment = " ") %>%
  summarise(
    N = n(),
    Mean = mean(area7_5k_nuc),
    Median = median(area7_5k_nuc),
    # Min = min(area7_5k_nuc),
    # Max = max(area7_5k_nuc),
    SD = sd(area7_5k_nuc),
    SEM = sd(area7_5k_nuc) / sqrt(length(area7_5k_nuc)),
    CV = sd(area7_5k_nuc) / mean(area7_5k_nuc) * 100
  ) %>% 
  gt() %>% 
  tab_header(title = "Analysis: Second Outermost (B+O) Edges") %>%
  fmt_number(decimals = 2)
