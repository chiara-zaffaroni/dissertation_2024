
# CM
# Mean and SD values for calculating Z' on selected thresholds
FF_Z_data_CM <- FF_plate %>% 
  filter(media_type %in% "CM") %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(lps, compound)) %>% 
  dplyr::select(Treatment,x17_5k_nuc, x20k_nuc, x25k_nuc) %>% 
  pivot_longer(cols = c(x17_5k_nuc, x20k_nuc, x25k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value)) %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value < 875)) %>%
  # filter(!(Treatment == "LPS.TAK3uM" & Threshold == "x2k_nuc" & Value > 1375)) %>%
  # filter(!(Treatment == "NoLPS.TAK3uM" & Threshold == "x2k_nuc" & Value < 900)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x3_5k_nuc" & Value > 1250)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x3_5k_nuc" & Value < 760)) %>%
  # filter(!(Treatment == "LPS.TAK3uM" & Threshold == "x3_5k_nuc" & Value > 1125)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x5k_nuc" & Value > 1125)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Mean = mean(Value),
            SD = sd(Value))


# Plot of all wrapped by threshold to check for outliers
FF_Z_CM_plot <- FF_Z_data_CM %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="CM")


# Group data by Threshold
grouped_data <- FF_Z_data_CM %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_mean, positive_sd, negative_mean, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_mean - positive_mean)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results_FF <- grouped_data %>%
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
      Mean[Treatment == "NoLPS.dmso"],
      SD[Treatment == "NoLPS.dmso"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "CM")


# MGM
# Mean and SD values for calculating Z' on selected thresholds
FF_Z_data_MGM <- FF_plate %>% 
  filter(media_type %in% "MGM") %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(lps, compound)) %>% 
  dplyr::select(Treatment, x2k_nuc, x3_5k_nuc, x5k_nuc) %>% 
  pivot_longer(cols = c(x2k_nuc, x3_5k_nuc, x5k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value)) %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value < 600)) %>%
  # filter(!(Treatment == "NoLPS.dmso" & Threshold == "x2k_nuc" & Value > 300)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x3_5k_nuc" & Value > 900)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x3_5k_nuc" & Value < 300)) %>%
  # filter(!(Treatment == "NoLPS.dmso" & Threshold == "x3_5k_nuc" & Value > 150)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x5k_nuc" & Value > 750)) %>%
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x5k_nuc" & Value < 225)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Mean = mean(Value),
            SD = sd(Value))


# Plot of all wrapped by threshold to check for outliers
FF_Z_MGM_plot <- FF_Z_data_MGM %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none") +
  labs(title ="MGM")


# Group data by Threshold
grouped_data_MGM <- FF_Z_data_MGM %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_mean, positive_sd, negative_mean, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_mean - positive_mean)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results_FF <- grouped_data_MGM %>%
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
      Mean[Treatment == "NoLPS.dmso"],
      SD[Treatment == "NoLPS.dmso"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "MGM")