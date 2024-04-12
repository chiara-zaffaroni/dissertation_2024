# 1. PREP
# load in packages
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

HH_plate_A3_original <- read_excel("hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx", sheet = 1, skip = 8) %>% 
  clean_names() %>% 
  # separate row and column
  separate(well_name, into = c("Row", "Column"), sep = -2) %>% 
  mutate(Column = as.numeric(Column)) %>% 
  arrange(Column) %>% 
  # add media type column
  mutate(Media_Type = "MGM") %>% 
  select(Media_Type, Row, Column, cell_nuclei_count_rod_h_mglia_dapi_cy5, everything()) %>% 
  rename("Nuclei_Count" ="cell_nuclei_count_rod_h_mglia_dapi_cy5") %>%
  select_if(~ !any(is.na(.))) %>%
  select(-c(plate_id, measurement_set_id, cell_object_id_rod_h_mglia_dapi_cy5)) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)

HH_LPS_map_A3 <- read_excel("hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

HH_platemap_A3 <- read_excel("hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx", sheet = 2, skip = 6) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(HH_LPS_map_A3$Value) %>% 
  rename(LPS_status = "...3")


HH_plate_A3 <- HH_plate_A3_original %>% 
  bind_cols(HH_platemap_A3[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(HH_plate_A3, "final_hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx")


HH_plate_A3 <- read_excel("final_hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx") %>% 
  clean_names()

HH_plate_A3_max <- HH_plate_A3 %>% 
  filter(compound == "dmso", lps_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

HH_A3_CV_analysis <- HH_plate_A3_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "Plate HH A3") %>%
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
# conclusion: use 1k or 2k threshold

# Z data prep (with edges)
LL_A1_Z_data <- LL_plate_A1 %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  filter(!Treatment %in% "noLPS.dmso") %>% 
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (with edges)
LL_A1_Z_plot <- LL_A1_Z_data %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30" = "3k Threshold", "cell_area_5k_rod_h_mglia_dapi_cy5_33" = "5k Threshold", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37" = "7.5k Threshold", "cell_area10k_rod_h_mglia_dapi_cy5_41" = "10k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title ="Plate LL A1 (with edges)", y = "CD38 Area/nuclei") 

# Group data by Threshold
grouped_data_LL_A1 <- LL_A1_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_LL_A1 <- grouped_data_LL_A1 %>%
  reframe(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "LPS.TAK3mM"],
      SD[Treatment == "LPS.TAK3mM"]
    ),
    Z_Prime_NoLPS_dmso = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "noLPS.dmso"],
      SD[Treatment == "noLPS.dmso"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "Plate LL A1 (with edges)")



# edgeless plate
HH_plate_A3_edgeless <- HH_plate_A3 %>% 
  filter(!column %in% c("1", "2", "23", "24")) %>% 
  filter(!row %in% c("A", "B", "P", "O")) 

HH_plate_A3_edgeless_max <- HH_plate_A3_edgeless %>% 
  filter(compound == "dmso", lps_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (edgeless)
HH_edgeless_CV_analysis <- HH_plate_A3_edgeless_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "Plate HH A3 (without edges)") %>%
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

# Z data prep (with edges)
HH_A3_Z_data <- HH_plate_A3 %>% 
  filter(compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5), names_to = "Threshold", values_to = "Value") %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  filter(!Treatment %in% "noLPS.TAK3uM") %>% 
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (no edges)
HH_A3_Z_plot <- HH_A3_Z_data %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot(colour = c("#00bb38", "#f88a82", "#72a5ff")) +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  scale_color_manual(values = c(noLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3uM='#72a5ff')) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5" = "1k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate HH A3 (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

comparisons_list <- list(c("noLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.TAK3uM"), c("noLPS.dmso", "LPS.TAK3uM"))
HH_A3_Z_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Group data by Threshold
grouped_data_HH_A3 <- HH_A3_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (no edges)
z_prime_results_HH_A3 <- grouped_data_HH_A3 %>%
  reframe(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "LPS.TAK3uM"],
      SD[Treatment == "LPS.TAK3uM"]
    ),
    Z_Prime_NoLPS_dmso = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "noLPS.dmso"],
      SD[Treatment == "noLPS.dmso"]
    )
  ) %>%  
  gt() %>% 
  tab_header(title = "Plate HH A3 (with edges)")


# Baricitinib 1K threshold
HH_A3_Z_data_Bari <- HH_plate_A3_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5, names_to = "Threshold", values_to = "Value") 

# Bari Plot
HH_A3_Z_Bari_plot <- HH_A3_Z_data_Bari %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Bari 0.3", "LPS.Bari 1", "LPS.Bari 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) +
  theme(legend.position = "none")+
  labs(title ="Plate HH A3 Baricitinib", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Bari 0.3"), c("LPS.dmso", "LPS.Bari 1"), c("LPS.dmso", "LPS.Bari 3"))
HH_A3_Z_Bari_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")

# Tacro 1K threshold
HH_A3_Z_data_Tacro <- HH_plate_A3_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5, names_to = "Threshold", values_to = "Value")

# Tacro Plot
HH_A3_Z_Tacro_plot <- HH_A3_Z_data_Tacro %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Tac 0.3", "LPS.Tac 1", "LPS.Tac 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate HH A3 Tacrolimus", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Tac 0.3"), c("LPS.dmso", "LPS.Tac 1"), c("LPS.dmso", "LPS.Tac 3"))
HH_A3_Z_Tacro_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Selinexor 1K threshold
HH_A3_Z_data_Sel <- HH_plate_A3_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Sel 3", "Sel 1", "Sel 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment, cell_area_1k_rod_h_mglia_dapi_cy5) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5, names_to = "Threshold", values_to = "Value")

# Selinexor Plot
HH_A3_Z_Sel_plot <- HH_A3_Z_data_Sel %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Sel 0.3", "LPS.Sel 1", "LPS.Sel 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate HH A3 Selinexor", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Sel 0.3"), c("LPS.dmso", "LPS.Sel 1"), c("LPS.dmso", "LPS.Sel 3"))
HH_A3_Z_Sel_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Everolimus 1K threshold
HH_A3_Z_data_Evero <- HH_plate_A3_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Ever 3", "Ever 1", "Ever 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5, names_to = "Threshold", values_to = "Value")

# Everolimus Plot
HH_A3_Z_Evero_plot <- HH_A3_Z_data_Evero %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Ever 0.3", "LPS.Ever 1", "LPS.Ever 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate HH A3 Everolimus", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Ever 0.3"), c("LPS.dmso", "LPS.Ever 1"), c("LPS.dmso", "LPS.Ever 3"))
HH_A3_Z_Evero_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


