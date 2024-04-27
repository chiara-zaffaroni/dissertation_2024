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

GG_plate_A1_original <- read_excel("hMglia GG A1 HighValue4 Widerange 10 threshold.txt.xlsx", sheet = 1, skip = 6) %>% 
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
          x15k_objects_total_rod_h_mglia_dapi_cy5_65, x3_5k_objects_total_rod_h_mglia_dapi_cy5_68, media_type)) %>% 
  select(Media_Type, Row, Column, nuclei_count_total_rod_h_mglia_dapi_cy5_14, everything()) %>% 
  rename("Nuclei_Count" ="nuclei_count_total_rod_h_mglia_dapi_cy5_14") %>%
  select(1:48) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)

GG_LPS_map_A1 <- read_excel("hMglia GG A1 HighValue4 Widerange 10 threshold.txt.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

GG_platemap_A1 <- read_excel("hMglia GG A1 HighValue4 Widerange 10 threshold.txt.xlsx", sheet = 2, skip = 2) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(GG_LPS_map_A1$Value) %>% 
  rename(LPS_status = "...3")


GG_plate_A1 <- GG_plate_A1_original %>% 
  bind_cols(GG_platemap_A1[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(GG_plate_A1, "final_hMglia GG A1 HighValue4 Widerange 10 threshold.xlsx")

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
GG_A1_max_area <- GG_plate_A1 %>%
  # filter for lps/compound combo first
  filter(LPS_status == "LPS" & Compound == "dmso") %>% 
  # group by media type
  group_by(Media_Type) %>% 
  select(1:5, contains("area"))

GG_A1_max_area_long <- GG_A1_max_area %>% 
  pivot_longer(cols = 6:16, names_to = "Threshold_using_area", values_to = "Value") 
# mutate(Value = as.numeric(Value))

GG_A1_max_area_long_plot <- GG_A1_max_area_long %>% 
  ggplot(aes(x = factor(Threshold_using_area, levels = c("cell_area_1k_rod_h_mglia_dapi_cy5_23", "cell_area_2k_rod_h_mglia_dapi_cy5_26", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30", "cell_area_5k_rod_h_mglia_dapi_cy5_33", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37", "cell_area10k_rod_h_mglia_dapi_cy5_41", "cell_area_12_5k_rod_h_mglia_dapi_cy5_44", "cell_area_15k_rod_h_mglia_dapi_cy5_48", "cell_area_17_5k_rod_h_mglia_dapi_cy5_51", "cell_area_20k_rod_h_mglia_dapi_cy5_55", "cell_area_25k_rod_h_mglia_dapi_cy5_59" )), y = Value, colour = Media_Type)) +  
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "hMglia GG Plate A1: Max (noLPS/dmso) signal using Area data",
       x = "Threshold",
       y = "CD38 Area/nuclei",
       colour = "Media Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("darkblue", "darkgreen"))

