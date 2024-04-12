II_plate_A2_original <- read_excel("hMglia II A2 MGM BariTacro Focus 48hrs.xlsx", sheet = 1, skip = 9) %>% 
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


II_LPS_map_A2 <- read_excel("hMglia II A2 MGM BariTacro Focus 48hrs.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

II_platemap_A2 <- read_excel("hMglia II A2 MGM BariTacro Focus 48hrs.xlsx", sheet = 2, skip = 3) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(II_LPS_map_A2$Value) %>% 
  rename(LPS_status = "...3")



II_plate_A2 <- II_plate_A2_original %>% 
  bind_cols(II_platemap_A2[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
# write_xlsx(II_plate_A2, "final_hMglia II A2 MGM BariTacro Focus 48hrs.xlsx")



# II_plate_A2_area <- II_plate_A2 %>% 
#   select(contains("area"))

II_plate_A2_area_long <- II_plate_A2 %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

II_plate_A2_area_plot <- II_plate_A2_area_long %>% 
  # edit naming of x axis tick labels and arrange in numerical order
  mutate(x_label = gsub("^cell_area10k.*", "cell_area10k", Variable)) %>%
  mutate(x_label = gsub("^(cell_area_[0-9]+(_[0-9]+)?k).*", "\\1", x_label)) %>%
  ggplot(aes(x = fct_inorder(x_label), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Plate II A2 (MGM 48h)",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")

II_plate_A2_max <- II_plate_A2 %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (with edges)
II_CV_analysis <- II_plate_A2_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "II A2 MGM 48h") %>%
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
# conclusion: use 1k or 2k or 3.5k threshold

# edgeless plate
II_plate_A2_edgeless <- II_plate_A2 %>% 
  filter(!Column %in% c("1", "2", "23", "24")) %>% 
  filter(!Row %in% c("A", "B", "P", "O")) 
II_plate_A2_edgeless_max <- II_plate_A2_edgeless %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (edgeless)
II_edgeless_CV_analysis <- II_plate_A2_edgeless_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "II A2 MGM 48h (without edges)") %>%
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
# conclusion: use 1k or 2k or 3.5k or 5k threshold

# Z data prep (with edges)
II_A2_Z_data <- II_plate_A2 %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30), names_to = "Threshold", values_to = "Value") %>% 
  # mutate(Value = as.numeric(Value)) %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (with edges)
II_A2_Z_plot <- II_A2_Z_data %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30" = "3k Threshold"))) +
  # theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate II A2 MGM 48h (with edges)", y = "CD38 Area/nuclei") 

# Group data by Threshold
grouped_data_II_A2 <- II_A2_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_II_A2 <- grouped_data_II_A2 %>%
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
  tab_header(title = "Plate II A2 MGM 48h (with edges)")



# Z data prep (no edges)
II_A2_Z_data_edgeless <- II_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30, cell_area_5k_rod_h_mglia_dapi_cy5_33) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26, cell_area_3_5k_rod_h_mglia_dapi_cy5_30, cell_area_5k_rod_h_mglia_dapi_cy5_33), names_to = "Threshold", values_to = "Value") %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (no edges)
II_A2_Z_plot_edgeless <- II_A2_Z_data_edgeless %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30" = "3k Threshold", "cell_area_5k_rod_h_mglia_dapi_cy5_33" = "5k Threshold"))) +
  # theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate II A2 MGM 48h (without edges)", y = "CD38 Area/nuclei") 

# Group data by Threshold
grouped_data_II_A2_edgeless <- II_A2_Z_data_edgeless %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (no edges)
z_prime_results_II_A2_edgeless <- grouped_data_II_A2_edgeless %>%
  reframe(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "LPS.TAK3mM"],
      SD[Treatment == "LPS.TAK3mM"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "Plate II A2 MGM 48h (without edges)")


# Baricitinib 1K threshold
II_A2_Z_data_Bari <- II_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Bari Plot
II_A2_Z_Bari_plot <- II_A2_Z_data_Bari %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate II A1 MGM 48h Baricitinib", y = "CD38 Area/nuclei")

# Tacro 1K threshold
II_A1_Z_data_Tacro <- II_plate_A1 %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Tacro Plot
II_A1_Z_Tacro_plot <- II_A1_Z_data_Tacro %>% 
  ggplot(aes(x=Treatment, y=Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate II A1 MGM 24h Tacrolimus", y = "CD38 Area/nuclei")
