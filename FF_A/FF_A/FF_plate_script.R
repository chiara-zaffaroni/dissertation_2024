library(readxl)
library(dplyr)
library(janitor)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(gt)
library(tseries)
library(tibble)
library(writexl)
library(stringr)
library(forcats)
library(ggsignif)

# read file
FF_plate_A1 <-  read_excel("hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx", sheet=1, skip=10) %>% 
  clean_names() %>% 
  select(-c(nuclei_count_total_rod_h_mglia_dapi_cy5, cell_area_2k_rod_h_mglia_dapi_cy5, cell_area_3_5k_rod_h_mglia_dapi_cy5, cell_area_5k_rod_h_mglia_dapi_cy5, cell_area_7_5k_rod_h_mglia_dapi_cy5, cell_area10k_rod_h_mglia_dapi_cy5, cell_area_12_5k_rod_h_mglia_dapi_cy5, cell_area_15k_rod_h_mglia_dapi_cy5, cell_area_17_5k_rod_h_mglia_dapi_cy5, cell_area_20k_rod_h_mglia_dapi_cy5, cell_area_25k_rod_h_mglia_dapi_cy5))

FF_plate_A1_CM <- FF_plate_A1 %>% 
  filter(!media_type %in% "MGM")

FF_plate_A1_MGM <- FF_plate_A1 %>% 
  filter(!media_type %in% "CM")

# CM: Area analysis
FF_plate_A1_CM_max <- FF_plate_A1_CM %>% 
  filter(compound == "dmso", lps == "LPS") %>% 
  pivot_longer(cols = contains("_nuc"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# MGM: Area analysis
FF_plate_A1_MGM_max <- FF_plate_A1_MGM %>% 
  filter(compound == "dmso", lps == "LPS") %>% 
  pivot_longer(cols = contains("_nuc"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CM: CV for max only 
FF_A1_CM_CV_analysis <- FF_plate_A1_CM_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "FF A1 CM 24h") %>%
  fmt_number(decimals = 2) %>% 
  tab_style(
    style = cell_text(
      size = "smaller",
      weight = "bold",
    ),
    locations = cells_body(columns = "CV")
  ) %>% 
  tab_style(
    locations = cells_body(
      columns = CV,
      rows = CV < 30
    ),
    style = list(cell_text(color = 'lightgreen')))
# conclusion: 

# MGM: CV for max only 
FF_A1_MGM_CV_analysis <- FF_plate_A1_MGM_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "FF A1 MGM 24h") %>%
  fmt_number(decimals = 2) %>% 
  tab_style(
    style = cell_text(
      size = "smaller",
      weight = "bold",
    ),
    locations = cells_body(columns = "CV")
  ) %>% 
  tab_style(
    locations = cells_body(
      columns = CV,
      rows = CV < 30
    ),
    style = list(cell_text(color = 'lightgreen')))
# conclusion: 


# CM: Z data prep (with edges)
FF_A1_CM_Z_data <- FF_plate_A1_CM %>% 
  filter(compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(lps, compound)) %>%
  dplyr::select(Treatment, x25k_nuc) %>% 
  pivot_longer(cols = c(x25k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  filter(!Treatment %in% "noLPS.dmso") %>% 
  summarise(Median = median(Value),
            SD = sd(Value)) 

grouped_data_FF_A1 <- FF_A1_CM_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_FF_A1 <- grouped_data_FF_A1 %>%
  reframe(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "LPS.TAK3uM"],
      SD[Treatment == "LPS.TAK3uM"]
    )
    #   Z_Prime_NoLPS_dmso = calculate_z_prime(
    #     Median[Treatment == "LPS.dmso"],
    #     SD[Treatment == "LPS.dmso"],
    #     Median[Treatment == "noLPS.dmso"],
    #     SD[Treatment == "noLPS.dmso"]
    #   )
  ) %>% 
  gt() %>% 
  tab_header(title = "Plate FF A1 CM 24h (with edges)")





# MEDIA EFFECT
FF_min_all_thresh <- FF_plate %>%
  # filter for lps/compound combo first
  filter(lps == "LPS" & compound == "TAK3uM") %>% 
  # group by media type
  group_by(media_type) %>% 
  select(media_type, lps, compound, x2k_nuc, x3_5k_nuc, x5k_nuc, x7_5k_nuc, x10k_nuc, x12_5k_nuc, x15k_nuc, x17_5k_nuc, x20k_nuc, x25k_nuc)

FF_min_all_thresh_long <- FF_min_all_thresh %>% 
  pivot_longer(cols = c(x2k_nuc, x3_5k_nuc, x5k_nuc, x7_5k_nuc, x10k_nuc, x12_5k_nuc, x15k_nuc, x17_5k_nuc, x20k_nuc, x25k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value))

FF_min_all_thresh_long_plot <- FF_min_all_thresh_long %>% 
  ggplot(aes(x = factor(Threshold, levels = c("x2k_nuc", "x3_5k_nuc", "x5k_nuc", "x7_5k_nuc", "x10k_nuc", "x12_5k_nuc", "x15k_nuc", "x17_5k_nuc", "x20k_nuc", "x25k_nuc")), y = Value, colour = media_type)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "hMglia FF Plate A3: Area all LPS/TAK (Min)",
       x = "Threshold",
       y = "CD38 Area/nuclei",
       colour = "Media Type") +
  theme_minimal() +
  scale_color_manual(values = c("darkblue", "darkgreen"))

# Define the order of thresholds
threshold_levels <- c("x2k_nuc", "x3_5k_nuc", "x5k_nuc", "x7_5k_nuc", "x10k_nuc", "x12_5k_nuc", "x15k_nuc", "x17_5k_nuc", "x20k_nuc", "x25k_nuc")

# JB Analysis
jarque_bera_results <- FF_max_all_thresh_long %>%
  mutate(Threshold = factor(Threshold, levels = threshold_levels)) %>%
  group_by(media_type, Threshold) %>%
  summarize(JB_Test_Statistic = jarque.bera.test(Value)$statistic,
            p_value = jarque.bera.test(Value)$p.value) %>%
  ungroup()

# Visualize the results
ggplot(jarque_bera_results, aes(x = Threshold, y = JB_Test_Statistic, color = media_type)) +
  geom_point() +
  geom_line() +
  labs(title = "Jarque-Bera Test Statistic for Each Threshold (by Media Type)",
       x = "Threshold",
       y = "Jarque-Bera Test Statistic",
       color = "Media Type")