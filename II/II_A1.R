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

II_plate_A1_original <- read_excel("hMglia II A1 MGM BariTacro Focus 24hrs.xlsx", sheet = 1, skip = 6) %>% 
  clean_names() %>% 
  # separate row and column
  separate(well_name, into = c("Row", "Column"), sep = -2) %>% 
  mutate(Column = as.numeric(Column)) %>% 
  arrange(Column) %>% 
  # add media type column
  mutate(Media_Type = "MGM") %>% 
  select(Media_Type, Row, Column, nuclei_count_total_rod_h_mglia_dapi_cy5_14, everything()) %>% 
  rename("Nuclei_Count" ="nuclei_count_total_rod_h_mglia_dapi_cy5_14") %>%
  # select_if(~ !any(is.na(.))) %>%
  select(-c(plate_id, measurement_set_id, cell_count_rod_h_mglia_dapi_cy5_4, cell_object_id_rod_h_mglia_dapi_cy5_5,
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
            x15k_objects_total_rod_h_mglia_dapi_cy5_65, x3_5k_objects_total_rod_h_mglia_dapi_cy5_68)) %>%
  select(1:48) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)


II_LPS_map_A1 <- read_excel("hMglia II A1 MGM BariTacro Focus 24hrs.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

II_platemap_A1 <- read_excel("hMglia II A1 MGM BariTacro Focus 24hrs.xlsx", sheet = 2, skip = 2) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(II_LPS_map_A1$Value) %>% 
  rename(LPS_status = "...3")



II_plate_A1 <- II_plate_A1_original %>% 
  bind_cols(II_platemap_A1[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
# write_xlsx(II_plate_A1, "final_hMglia II A1 MGM BariTacro Focus 24hrs.xlsx")

II_plate_A1 <- read_excel("final_hMglia II A1 MGM BariTacro Focus 24hrs.xlsx") %>% 
  clean_names()

II_plate_A1_area <- II_plate_A1 %>% 
  select(contains("area"))
  
II_plate_A1_area_long <- II_plate_A1 %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

II_plate_A1_area_plot <- II_plate_A1_area_long %>% 
  # edit naming of x axis tick labels and arrange in numerical order
  mutate(x_label = gsub("^cell_area10k.*", "cell_area10k", Variable)) %>%
  mutate(x_label = gsub("^(cell_area_[0-9]+(_[0-9]+)?k).*", "\\1", x_label)) %>%  
  ggplot(aes(x = fct_inorder(x_label), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Plate II A1 (MGM 24h)",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")

II_plate_A1_max <- II_plate_A1 %>% 
  filter(compound == "dmso", lps_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only
II_CV_analysis <- II_plate_A1_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "II A1 MGM 24h") %>%
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

# edgeless data
II_plate_A1_edgeless <- II_plate_A1 %>% 
  filter(!column %in% c("1", "2", "23", "24")) %>% 
  filter(!row %in% c("A", "B", "P", "O")) 

II_plate_A1_edgeless_max <- II_plate_A1_edgeless %>% 
  filter(compound == "dmso", lps_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)  

II_edgeless_CV_analysis <- II_plate_A1_edgeless_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "II A1 MGM 24h (no edges)") %>%
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




II_A1_Z_data <- II_plate_A1 %>% 
  filter(compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23, ) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23), names_to = "Threshold", values_to = "Value") %>% 
  # mutate(Value = as.numeric(Value)) %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  filter(!Treatment %in% "noLPS.TAK3uM") %>% 
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot
II_A1_Z_plot <- II_A1_Z_data %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot(colour = c("#00bb38", "#f88a82", "#72a5ff")) +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  scale_color_manual(values = c(noLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3mM='#72a5ff')) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate II A1 MGM 24h (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

comparisons_list <- list(c("noLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.TAK3mM"), c("noLPS.dmso", "LPS.TAK3mM"))
II_A1_Z_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")



# Group data by Threshold
grouped_data_II_A1 <- II_A1_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations
z_prime_results_II_A1 <- grouped_data_II_A1 %>%
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
  tab_header(title = "Plate II A1 MGM 24h")


# Baricitinib 1K threshold
II_A1_Z_data_Bari <- II_plate_A1_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Bari Plot
II_A1_Z_Bari_plot <- II_A1_Z_data_Bari %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Bari 0.3", "LPS.Bari 1", "LPS.Bari 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate II A1 MGM 24h Baricitinib", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Bari 0.3"), c("LPS.dmso", "LPS.Bari 1"), c("LPS.dmso", "LPS.Bari 3"))
II_A1_Z_Bari_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")

# Tacro 1K threshold
II_A1_Z_data_Tacro <- II_plate_A1_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Tacro Plot
II_A1_Z_Tacro_plot <- II_A1_Z_data_Tacro %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Tac 0.3", "LPS.Tac 1", "LPS.Tac 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate II A1 MGM 24h Tacrolimus", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Tac 0.3"), c("LPS.dmso", "LPS.Tac 1"), c("LPS.dmso", "LPS.Tac 3"))
II_A1_Z_Tacro_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Selinexor 1K threshold
II_A1_Z_data_Sel <- II_plate_A1_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Sel 3", "Sel 1", "Sel 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment, cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Selinexor Plot
II_A1_Z_Sel_plot <- II_A1_Z_data_Sel %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Sel 0.3", "LPS.Sel 1", "LPS.Sel 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate II A1 Selinexor", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Sel 0.3"), c("LPS.dmso", "LPS.Sel 1"), c("LPS.dmso", "LPS.Sel 3"))
II_A1_Z_Sel_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Everolimus 1K threshold
II_A1_Z_data_Evero <- II_plate_A1_edgeless %>% 
  filter(compound %in% c("dmso", "TAK3mM", "Ever 3", "Ever 1", "Ever 0.3")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Everolimus Plot
II_A1_Z_Evero_plot <- II_A1_Z_data_Evero %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Ever 0.3", "LPS.Ever 1", "LPS.Ever 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate II A1 Everolimus", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Ever 0.3"), c("LPS.dmso", "LPS.Ever 1"), c("LPS.dmso", "LPS.Ever 3"))
II_A1_Z_Evero_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# _________________________
# do nuclei count graph for each drug

II_A1_nuclei_Tacro <- II_plate_A1 %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,Nuclei_Count)

nuclei_Tacro_plot <- II_A1_nuclei_Tacro %>% 
  ggplot(aes(x=Treatment, y=Nuclei_Count, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  theme(legend.position = "none")+
  labs(title ="Tacrolimus Nuclei Counts", y = "Nuclei count")
