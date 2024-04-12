LL_plate_A2_original <- read_excel("hMglia LL A2 IL1B TacroBari round4 15.3.24.xlsx", sheet = 1, skip = 8) %>% 
  clean_names() %>% 
  separate(well_name, into = c("row", "column"), sep = -2) %>% 
  mutate(column = as.numeric(column)) %>% 
  arrange(column) %>% 
  select(row, column, nuclei_count_total_rod_h_mglia_dapi_cy5_14, everything()) %>% 
  rename("Nuclei_Count" ="nuclei_count_total_rod_h_mglia_dapi_cy5_14") %>%
  select(-c(cell_10k_objects_rod_h_mglia_dapi_cy5_13, cell_10k_objects_rod_h_mglia_dapi_cy5_72,cell_12_5k_objects_rod_h_mglia_dapi_cy5_106,
            cell_12_5k_objects_rod_h_mglia_dapi_cy5_47, cell_15k_objects_rod_h_mglia_dapi_cy5_66,
            cell_15k_objects_rod_h_mglia_dapi_cy5_7, cell_17_5k_objects_rod_h_mglia_dapi_cy5_113,
            cell_17_5k_objects_rod_h_mglia_dapi_cy5_54, cell_1k_objects_rod_h_mglia_dapi_cy5_11, 
            cell_1k_objects_rod_h_mglia_dapi_cy5_70, cell_20k_objects_rod_h_mglia_dapi_cy5_117,
            cell_20k_objects_rod_h_mglia_dapi_cy5_58, cell_25k_objects_rod_h_mglia_dapi_cy5_121,
            cell_25k_objects_rod_h_mglia_dapi_cy5_62, cell_2k_objects_rod_h_mglia_dapi_cy5_29,
            cell_2k_objects_rod_h_mglia_dapi_cy5_88, cell_3_5k_objects_rod_h_mglia_dapi_cy5_12,
            cell_3_5k_objects_rod_h_mglia_dapi_cy5_71, cell_5k_objects_rod_h_mglia_dapi_cy5_36,
            cell_5k_objects_rod_h_mglia_dapi_cy5_95, cell_7_5k_objects_rod_h_mglia_dapi_cy5_40,
            cell_7_5k_objects_rod_h_mglia_dapi_cy5_99,cell_nuclei_count_rod_h_mglia_dapi_cy5_22,
            x15k_objects_total_rod_h_mglia_dapi_cy5_65, x3_5k_objects_total_rod_h_mglia_dapi_cy5_68, plate_id, measurement_set_id, cell_object_id_rod_h_mglia_dapi_cy5_5, cell_count_rod_h_mglia_dapi_cy5_4)) %>%
  select(1:47) %>% 
  mutate_at(vars(4:47), ~ . / Nuclei_Count)

LL_LPS_map_A2 <- read_excel("hMglia LL A2 IL1B TacroBari round4 15.3.24.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

LL_platemap_A2 <- read_excel("hMglia LL A2 IL1B TacroBari round4 15.3.24.xlsx", sheet = 2) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(LL_LPS_map_A2$Value) %>% 
  rename(LPS_status = "...3")

LL_plate_A2 <- LL_plate_A2_original %>% 
  bind_cols(LL_platemap_A2[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(row, column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(LL_plate_A2, "final_hMglia LL A2 IL1B TacroBari round4 15.3.24.xlsx")     


# Area analysis

# LL_plate_A1_area <- LL_plate_A1 %>% 
#   select(contains("area"))
# 
# LL_plate_A1_area_long <- LL_plate_A1 %>% 
#   pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
#   select(Variable, Value)

LL_plate_A2_max <- LL_plate_A2 %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

LL_plate_A2_area_plot <- LL_plate_A2_max %>% 
  # edit naming of x axis tick labels and arrange in numerical order
  mutate(x_label = gsub("^cell_area10k.*", "cell_area10k", Variable)) %>%
  mutate(x_label = gsub("^(cell_area_[0-9]+(_[0-9]+)?k).*", "\\1", x_label)) %>%  
  ggplot(aes(x = fct_inorder(x_label), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Plate LL A2",
       x = "Threshold",
       y = "IL1B Area/nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")


# CV for max only (with edges)
LL_A2_CV_analysis <- LL_plate_A2_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "Plate LL A2") %>%
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
LL_A2_Z_data <- LL_plate_A2 %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26, ) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (with edges)
LL_A2_Z_plot <- LL_A2_Z_data %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold"))) +
  # theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title ="Plate LL A2 (with edges)", y = "IL1B Area/nuclei") 

# Group data by Threshold
grouped_data_LL_A2 <- LL_A2_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_LL_A2 <- grouped_data_LL_A2 %>%
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
  tab_header(title = "Plate LL A2 (with edges)")



# edgeless plate
LL_plate_A2_edgeless <- LL_plate_A2 %>% 
  filter(!column %in% c("1", "2", "23", "24")) %>% 
  filter(!row %in% c("A", "B", "P", "O")) 
LL_plate_A2_edgeless_max <- LL_plate_A2_edgeless %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (edgeless)
LL_edgeless_CV_analysis <- LL_plate_A2_edgeless_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "Plate LL A2 (without edges)") %>%
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

# Z data prep (no edges)
LL_A2_Z_data_edgeless <- LL_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23, cell_area_2k_rod_h_mglia_dapi_cy5_26), names_to = "Threshold", values_to = "Value") %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (no edges)
LL_A2_Z_plot_edgeless <- LL_A2_Z_data_edgeless %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate LL A2 (without edges)", y = "IL1B Area/nuclei") 

# Group data by Threshold
grouped_data_LL_A2_edgeless <- LL_A2_Z_data_edgeless %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (no edges)
z_prime_results_LL_A2_edgeless <- grouped_data_LL_A2_edgeless %>%
  reframe(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "LPS.TAK3mM"],
      SD[Treatment == "LPS.TAK3mM"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "Plate LL A2 (without edges)")


# Baricitinib 1K threshold
LL_A2_Z_data_Bari <- LL_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value") 

# Bari Plot
LL_A2_Z_Bari_plot <- LL_A2_Z_data_Bari %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Bari 0.3", "LPS.Bari 1", "LPS.Bari 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) +
  theme(legend.position = "none")+
  labs(title ="Plate LL A2 Baricitinib", y = "IL1B Area/nuclei")

# checking assumptions for 1 factor GLM
# Specify reference level for the category variable
LL_A2_Z_data_Bari$Treatment <- relevel(LL_A2_Z_data_Bari$Treatment, ref = "LPS.dmso")
Bari_lm <- lm(Value ~ factor(Treatment), data = LL_A2_Z_data_Bari)
plot(Bari_lm)
summary(Bari_lm)
anova(Bari_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Bari 0.3"), c("LPS.dmso", "LPS.Bari 1"), c("LPS.dmso", "LPS.Bari 3"))
LL_A2_Z_Bari_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# # calculating the mean values for each type of Treatment
# Bari_emmeans <- emmeans(Bari_lm, "Treatment")
# Bari_pairs <- pairs(Bari_emmeans)
# confint(Bari_pairs)

# Tacro 1K threshold
LL_A2_Z_data_Tacro <- LL_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Tacro Plot
LL_A2_Z_Tacro_plot <- LL_A2_Z_data_Tacro %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Tac 0.3", "LPS.Tac 1", "LPS.Tac 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate LL A2 Tacrolimus", y = "IL1B Area/nuclei")

LL_A2_Z_data_Tacro$Treatment <- relevel(LL_A2_Z_data_Tacro$Treatment, ref = "LPS.dmso")
Tacro_lm <- lm(Value ~ factor(Treatment), data = LL_A2_Z_data_Tacro)
plot(Tacro_lm)
summary(Tacro_lm)
anova(Tacro_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Tac 0.3"), c("LPS.dmso", "LPS.Tac 1"), c("LPS.dmso", "LPS.Tac 3"))
LL_A2_Z_Tacro_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Selinexor 1K threshold
LL_A2_Z_data_Sel <- LL_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Sel 3", "Sel 1", "Sel 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Selinexor Plot
LL_A2_Z_Sel_plot <- LL_A2_Z_data_Sel %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Sel 0.3", "LPS.Sel 1", "LPS.Sel 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate LL A2 Selinexor", y = "IL1B Area/nuclei")

Sel_lm <- lm(Value ~ factor(Treatment), data = LL_A2_Z_data_Sel)
plot(Sel_lm)
summary(Sel_lm)
anova(Sel_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Sel 0.3"), c("LPS.dmso", "LPS.Sel 1"), c("LPS.dmso", "LPS.Sel 3"))
LL_A2_Z_Sel_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Everolimus 1K threshold
LL_A2_Z_data_Evero <- LL_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Ever 3", "Ever 1", "Ever 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Everolimus Plot
LL_A2_Z_Evero_plot <- LL_A2_Z_data_Evero %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Ever 0.3", "LPS.Ever 1", "LPS.Ever 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate LL A2 Everolimus", y = "IL1B Area/nuclei")

Evero_lm <- lm(Value ~ factor(Treatment), data = LL_A2_Z_data_Evero)
plot(Evero_lm)
summary(Evero_lm)
anova(Evero_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Ever 0.3"), c("LPS.dmso", "LPS.Ever 1"), c("LPS.dmso", "LPS.Ever 3"))
LL_A2_Z_Evero_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")