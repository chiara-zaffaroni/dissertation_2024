# Mean and SD values for calculating Z'
EE_Z_data_edgefree <- EE_plate_edgefree %>% 
  filter(media_type %in% "CM") %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(lps, compound)) %>% 
  filter(!Treatment %in% "NoLPS.TAK3uM") %>% 
  dplyr::select(Treatment, area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc) %>% 
  pivot_longer(cols = c(area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value)) %>%   
  # remove outliers
  filter(!(Threshold == "area5k_nuc" & Value > 600)) %>%
  filter(!(Treatment == "LPS.dmso" & Threshold == "area5k_nuc" & Value < 100)) %>%
  filter(!(Treatment == "NoLPS.dmso" & Threshold == "area5k_nuc" & Value > 125)) %>%
  filter(!(Threshold == "area7_5k_nuc" & Value < 10)) %>%
  filter(!(Threshold == "area7_5k_nuc" & Value > 300)) %>%
  filter(!(Treatment == "LPS.dmso" & Threshold == "area10k_nuc" & Value < 50)) %>%
  filter(!(Treatment == "LPS.TAK3uM" & Threshold == "area10k_nuc" & Value > 150)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Mean = mean(Value),
            SD = sd(Value))
  # gt()

# Plot of all wrapped by threshold to check for outliers
plot <- EE_Z_data_edgefree %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")
  

# Group data by Threshold
grouped_data <- EE_Z_data_edgefree %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_mean, positive_sd, negative_mean, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_mean - positive_mean)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results <- grouped_data %>%
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
  gt()

# z_prime_values <- EE_Z_data_edgefree %>%
#   group_by(Threshold) %>%
#   summarise(
#     Positive_Control_Mean = Mean[Treatment == "LPS.dmso"],
#     Positive_Control_SD = SD[Treatment == "LPS.dmso"],
#     Negative_Control_Mean = Mean[Treatment %in% c("LPS.TAK3uM", "NoLPS.TAK3uM")],
#     Negative_Control_SD = SD[Treatment %in% c("LPS.TAK3uM", "NoLPS.TAK3uM")],
#     Z_Prime = 1 - (3 * (Negative_Control_SD + Positive_Control_SD)) / abs(Negative_Control_Mean - Positive_Control_Mean)
#   ) %>%
#   filter(!is.infinite(Z_Prime) & !is.nan(Z_Prime))


