# CM: want to select 1K, 2K thresholds for area
GG_Z_data_CM <- GG_plate_A1 %>% 
  filter(Media_Type == "CM") %>% 
  filter(Compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, cell_area_1k_rod_h_mglia_dapi_cy5_23,  cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30, cell_area_5k_rod_h_mglia_dapi_cy5_33) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23,  cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30, cell_area_5k_rod_h_mglia_dapi_cy5_33), names_to = "Threshold", values_to = "Value") %>% 
  # mutate(Value = as.numeric(Value)) %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Mean = mean(Value),
            SD = sd(Value))
  
# Plot of all wrapped by threshold to check for outliers
GG_Z_data_CM_plot <- GG_Z_data_CM %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="CM")

# Group data by Threshold
grouped_data <- GG_Z_data_CM %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_mean, positive_sd, negative_mean, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_mean - positive_mean)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results_GG <- grouped_data %>%
  summarize(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Mean[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Mean[Treatment == "LPS.TAK3uM"],
      SD[Treatment == "LPS.TAK3uM"]
    ),
    Z_Prime_NoLPS_dmso = calculate_z_prime(
      Mean[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Mean[Treatment == "noLPS.dmso"],
      SD[Treatment == "noLPS.dmso"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "CM")



# MGM: want to select 1K, 2K, 3.5K, 5K thresholds for area
GG_Z_data_MGM <- GG_plate_A1 %>% 
  filter(Media_Type %in% "MGM") %>% 
  filter(Compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, cell_area_1k_rod_h_mglia_dapi_cy5_23,  cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23,  cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30), names_to = "Threshold", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value)) %>%
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Mean = mean(Value),
            SD = sd(Value))

# Plot of all wrapped by threshold to check for outliers
GG_Z_data_MGM_plot <- GG_Z_data_MGM %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="MGM", y = "CD38 Area/nuclei")

# Group data by Threshold
grouped_data_MGM <- GG_Z_data_MGM %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_mean, positive_sd, negative_mean, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_mean - positive_mean)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results_GG <- GG_Z_data_MGM %>%
  summarize(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Mean[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Mean[Treatment == "LPS.TAK3uM"],
      SD[Treatment == "LPS.TAK3uM"]
    ),
    Z_Prime_NoLPS_dmso = calculate_z_prime(
      Mean[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Mean[Treatment == "noLPS.dmso"],
      SD[Treatment == "noLPS.dmso"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "MGM")
