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

GG_plate_A3_original <- read_excel("hMglia GG A3 MinMax LPS dmso vs TAK242.xlsx", sheet = 1, skip = 9) %>% 
  clean_names() %>% 
  # separate row and column
  separate(well_name, into = c("Row", "Column"), sep = -2) %>% 
  mutate(Column = as.numeric(Column)) %>% 
  arrange(Column) %>% 
  # add media type column
  mutate(Media_Type = case_when(
    row_number() <= 192 ~ "CM",
    row_number() > 192 ~ "MGM"
  )) %>% 
  select(-c(measurement_set_id, cell_count_rod_h_mglia_dapi_cy5_4, cell_object_id_rod_h_mglia_dapi_cy5_5,
            cell_10k_objects_rod_h_mglia_dapi_cy5_13, cell_10k_objects_rod_h_mglia_dapi_cy5_72,cell_12_5k_objects_rod_h_mglia_dapi_cy5_106,
            cell_12_5k_objects_rod_h_mglia_dapi_cy5_47, cell_15k_objects_rod_h_mglia_dapi_cy5_66,
            cell_15k_objects_rod_h_mglia_dapi_cy5_7, cell_17_5k_objects_rod_h_mglia_dapi_cy5_113,
            cell_17_5k_objects_rod_h_mglia_dapi_cy5_54, cell_1k_objects_rod_h_mglia_dapi_cy5_11, 
            cell_1k_objects_rod_h_mglia_dapi_cy5_70, cell_20k_objects_rod_h_mglia_dapi_cy5_117,
            cell_20k_objects_rod_h_mglia_dapi_cy5_58, cell_25k_objects_rod_h_mglia_dapi_cy5_121,
            cell_25k_objects_rod_h_mglia_dapi_cy5_62, cell_2k_objects_rod_h_mglia_dapi_cy5_29,
            cell_2k_objects_rod_h_mglia_dapi_cy5_88, cell_3_5k_objects_rod_h_mglia_dapi_cy5_12,
            cell_3_5k_objects_rod_h_mglia_dapi_cy5_71, cell_5k_objects_rod_h_mglia_dapi_cy5_36,
            cell_5k_objects_rod_h_mglia_dapi_cy5_95, cell_7_5k_objects_rod_h_mglia_dapi_cy5_40,
            cell_7_5k_objects_rod_h_mglia_dapi_cy5_99, cell_count_rod_h_mglia_dapi_cy5_63,
            cell_nuclei_count_rod_h_mglia_dapi_cy5_22, cell_object_id_rod_h_mglia_dapi_cy5_64,
            x15k_objects_total_rod_h_mglia_dapi_cy5_65, x3_5k_objects_total_rod_h_mglia_dapi_cy5_68, plate_id)) %>% 
  select(Media_Type, Row, Column, nuclei_count_total_rod_h_mglia_dapi_cy5_14, everything()) %>% 
  rename("Nuclei_Count" ="nuclei_count_total_rod_h_mglia_dapi_cy5_14") %>%
  select(1:48) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)



GG_LPS_map_A3 <- read_excel("hMglia GG A3 MinMax LPS dmso vs TAK242.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

GG_platemap_A3 <- read_excel("hMglia GG A3 MinMax LPS dmso vs TAK242.xlsx", sheet = 2, skip = 3) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(GG_LPS_map_A3$Value) %>% 
  rename(LPS_status = "...3")



GG_plate_A3 <- GG_plate_A3_original %>% 
  bind_cols(GG_platemap_A1[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(GG_plate_A3, "final_hMglia GG A3 MinMax LPS dmso vs TAK242.xlsx")

# _______



GG_plate_A3_CM <- read_excel("final_hMglia GG A3 MinMax LPS dmso vs TAK242.xlsx", sheet = 3) %>% 
  clean_names()

GG_plate_A3_CM_max <- GG_plate_A3_CM %>% 
  filter(compound == "dmso", lps_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (with edges)
GG_plate_A3_CM_CV_analysis <- GG_plate_A3_CM_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "GG A3 CM 24h") %>%
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
# conclusion: 1k

# Z data prep (with edges)
GG_A3_CM_Z_data <- GG_plate_A3_CM %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>%
  dplyr::select(Treatment, area_2k) %>% 
  pivot_longer(cols = area_2k, names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (with edges)
GG_A3_CM_Z_plot <- GG_A3_CM_Z_data %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("area_2k" = "2K Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate GG A3 CM 24h (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(100, 1600))

# Group data by Threshold
grouped_data_GG_CM_A3 <- GG_A3_CM_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_GG_A3 <- grouped_data_GG_CM_A3 %>%
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
  tab_header(title = "Plate GG A3 CM 24h (with edges)")


GG_plate_A3_MGM <- read_excel("final_hMglia GG A3 MinMax LPS dmso vs TAK242.xlsx", sheet = 4) %>% 
  clean_names()

GG_plate_A3_MGM_max <- GG_plate_A3_MGM %>% 
  filter(compound == "dmso", lps_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (with edges)
GG_plate_A3_MGM_CV_analysis <- GG_plate_A3_MGM_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "GG A3 MGM 24h") %>%
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
# conclusion: 1k

# Z data prep (with edges)
GG_A3_MGM_Z_data <- GG_plate_A3_MGM %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>%
  dplyr::select(Treatment, area_2k) %>% 
  pivot_longer(cols = area_2k, names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (with edges)
GG_A3_MGM_Z_plot <- GG_A3_MGM_Z_data %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("area_2k" = "2K Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate GG A3 MGM 24h (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(100, 1600))

# Group data by Threshold
grouped_data_GG_MGM_A3 <- GG_A3_MGM_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_GG_MGM_A3 <- grouped_data_GG_MGM_A3 %>%
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
  tab_header(title = "Plate GG A3 MGM 24h (with edges)")





# MIN AREA SIGNAL FOR BOTH MEDIA
GG_A1_min_area <- GG_plate_A1 %>%
  # filter for lps/compound combo first
  filter(LPS_status == "LPS" & Compound == "TAK3uM") %>% 
  # group by media type
  group_by(Media_Type) %>% 
  select(1:5, contains("area"))

GG_A1_min_area_long <- GG_A1_min_area %>% 
  pivot_longer(cols = 6:16, names_to = "Threshold_using_area", values_to = "Value") 
# mutate(Value = as.numeric(Value))

GG_A1_min_area_long_plot <- GG_A1_min_area_long %>% 
  ggplot(aes(x = factor(Threshold_using_area, levels = c("cell_area_1k_rod_h_mglia_dapi_cy5_23", "cell_area_2k_rod_h_mglia_dapi_cy5_26", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30", "cell_area_5k_rod_h_mglia_dapi_cy5_33", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37", "cell_area10k_rod_h_mglia_dapi_cy5_41", "cell_area_12_5k_rod_h_mglia_dapi_cy5_44", "cell_area_15k_rod_h_mglia_dapi_cy5_48", "cell_area_17_5k_rod_h_mglia_dapi_cy5_51", "cell_area_20k_rod_h_mglia_dapi_cy5_55", "cell_area_25k_rod_h_mglia_dapi_cy5_59" )), y = Value, colour = Media_Type)) +  
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "hMglia GG Plate A1: Min (LPS/TAK) signal using Area data",
       x = "Threshold",
       y = "CD38 Area/nuclei",
       colour = "Media Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  # theme_minimal() +
  scale_color_manual(values = c("darkblue", "darkgreen"))


# MAX AREA SIGNAL FOR BOTH MEDIA
GG_A3_max_area <- GG_plate_A3 %>%
  # filter for lps/compound combo first
  filter(LPS_status == "LPS" & Compound == "dmso") %>% 
  # group by media type
  group_by(Media_Type) %>% 
  select(1:5, contains("area"))

GG_A3_max_area_long <- GG_A3_max_area %>% 
  pivot_longer(cols = 6:16, names_to = "Threshold_using_area", values_to = "Value") 
# mutate(Value = as.numeric(Value))

GG_A3_max_area_long_plot <- GG_A3_max_area_long %>% 
  ggplot(aes(x = factor(Threshold_using_area, levels = c("cell_area_1k_rod_h_mglia_dapi_cy5_23", "cell_area_2k_rod_h_mglia_dapi_cy5_26", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30", "cell_area_5k_rod_h_mglia_dapi_cy5_33", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37", "cell_area10k_rod_h_mglia_dapi_cy5_41", "cell_area_12_5k_rod_h_mglia_dapi_cy5_44", "cell_area_15k_rod_h_mglia_dapi_cy5_48", "cell_area_17_5k_rod_h_mglia_dapi_cy5_51", "cell_area_20k_rod_h_mglia_dapi_cy5_55", "cell_area_25k_rod_h_mglia_dapi_cy5_59" )), y = Value, colour = Media_Type)) +  
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "hMglia GG Plate A3: Max (LPS/dmso) signal using Area data",
       x = "Threshold",
       y = "CD38 Area/nuclei",
       colour = "Media Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("darkblue", "darkgreen"))


# MIN AREA SIGNAL FOR BOTH MEDIA
GG_A3_min_area <- GG_plate_A3 %>%
  # filter for lps/compound combo first
  filter(LPS_status == "LPS" & Compound == "TAK3uM") %>% 
  # group by media type
  group_by(Media_Type) %>% 
  select(1:5, contains("area"))

GG_A3_min_area_long <- GG_A3_min_area %>% 
  pivot_longer(cols = 6:16, names_to = "Threshold_using_area", values_to = "Value") 
# mutate(Value = as.numeric(Value))

GG_A3_min_area_long_plot <- GG_A3_min_area_long %>% 
  ggplot(aes(x = factor(Threshold_using_area, levels = c("cell_area_1k_rod_h_mglia_dapi_cy5_23", "cell_area_2k_rod_h_mglia_dapi_cy5_26", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30", "cell_area_5k_rod_h_mglia_dapi_cy5_33", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37", "cell_area10k_rod_h_mglia_dapi_cy5_41", "cell_area_12_5k_rod_h_mglia_dapi_cy5_44", "cell_area_15k_rod_h_mglia_dapi_cy5_48", "cell_area_17_5k_rod_h_mglia_dapi_cy5_51", "cell_area_20k_rod_h_mglia_dapi_cy5_55", "cell_area_25k_rod_h_mglia_dapi_cy5_59" )), y = Value, colour = Media_Type)) +  
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "hMglia GG Plate A3: Min (LPS/TAK-242) signal using Area data",
       x = "Threshold",
       y = "CD38 Area/nuclei",
       colour = "Media Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  # theme_minimal() +
  scale_color_manual(values = c("darkblue", "darkgreen"))


# CVs for max MGM plate 
GG_CV_analysis_MGM_A3 <- GG_A3_max_area_long %>% 
  filter(Media_Type == "MGM") %>% 
  group_by(Threshold_using_area) %>% 
  select(Threshold_using_area, Value) %>% 
  # mutate(Mean = mean(Value))
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "MGM media") %>%
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
    style = list(cell_text(color = 'lightgreen'))
  )
# conclusion: 


# MGM: want to select 1K, 2K, 3.5K, 5K thresholds for area
GG_Z_data_MGM_A3 <- GG_plate_A3 %>% 
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
GG_Z_data_MGM_A3_plot <- GG_Z_data_MGM_A3 %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="MGM", y = "CD38 Area/nuclei")

# # Group data by Threshold
# grouped_data_MGM_A3 <- GG_Z_data_MGM_A3 %>%
#   group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_mean, positive_sd, negative_mean, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_mean - positive_mean)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results_GG_A3 <- GG_Z_data_MGM_A3 %>%
  summarize(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Mean[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Mean[Treatment == "LPS.TAK3uM"],
      SD[Treatment == "LPS.TAK3uM"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "MGM")
