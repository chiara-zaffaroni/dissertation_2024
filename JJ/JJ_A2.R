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

JJ_plate_A2_original <- read_excel("hMglia JJ A2.xlsx", sheet = 1) %>% 
  clean_names() %>% 
  arrange(column, row) %>%
  select(media_type, row, column, nuclei_count_total_rod_h_mglia_dapi_cy5_14, everything()) %>% 
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
            x15k_objects_total_rod_h_mglia_dapi_cy5_65, x3_5k_objects_total_rod_h_mglia_dapi_cy5_68, compound, lps_status)) %>%
  select(1:48) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)


JJ_LPS_map_A2 <- read_excel("hMglia JJ A2.xlsx", sheet = 4) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

JJ_platemap_A2 <- read_excel("hMglia JJ A2.xlsx", sheet = 3, skip = 5) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(JJ_LPS_map_A2$Value) %>% 
  rename(LPS_status = "...3")

JJ_plate_A2 <- JJ_plate_A2_original %>% 
  bind_cols(JJ_platemap_A2[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(media_type, row, column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
# write_xlsx(JJ_plate_A2, "final_hMglia JJ A2.xlsx")     


# Area analysis
JJ_plate_A2_max <- JJ_plate_A2 %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

JJ_plate_A2_area_plot <- JJ_plate_A2_max %>% 
  # edit naming of x axis tick labels and arrange in numerical order
  mutate(x_label = gsub("^cell_area10k.*", "cell_area10k", Variable)) %>%
  mutate(x_label = gsub("^(cell_area_[0-9]+(_[0-9]+)?k).*", "\\1", x_label)) %>%  
  ggplot(aes(x = fct_inorder(x_label), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Plate JJ A2 (MGM 24h)",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")



# CV for max only (with edges)
JJ_A2_CV_analysis <- JJ_plate_A2_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "JJ A2 MGM 24h") %>%
  fmt_number(decimals = 2) %>% 
  tab_style(
    style = cell_text(
      size = "smaller",
      weight = "bold",
    ),
    locations = cells_body(columns = "CV")
  ) %>% 
  # tab_style(
  #   style = cell_fill(color = "lightgreen"),
  #   locations = cells_body(columns = "CV")
  # ) %>% 
  tab_style(
    locations = cells_body(
      columns = CV,
      rows = CV < 30
    ),
    style = list(cell_text(color = 'lightgreen')))
# conclusion: use 1k threshold

# Z data prep (with edges)
JJ_A2_Z_data <- JJ_plate_A2 %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = c(cell_area_1k_rod_h_mglia_dapi_cy5_23), names_to = "Threshold", values_to = "Value") %>% 
  filter(!Treatment %in% "noLPS.dmso") %>% 
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (with edges)
JJ_A2_Z_plot <- JJ_A2_Z_data %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate JJ A2 CM 24h (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 1000))

# JJ_A1_Z_data <- mutate(JJ_A1_Z_data, dataset = "Dataset 1")
# JJ_A2_Z_data <- mutate(JJ_A2_Z_data, dataset = "Dataset 2")
# combined_data <- bind_rows(JJ_A1_Z_data, JJ_A2_Z_data)
# 
# ggplot(combined_data, aes(x = Treatment, y = Value)) +
#   geom_boxplot() +
#   labs(title = "Comparison of Phenotypes",
#        x = "Cell line", y = "% Number of individual cells") +
#   facet_wrap(~dataset, scales = "free", ncol = 2) +
#   theme(axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20))
                                    

# Group data by Threshold
grouped_data_JJ_A2 <- JJ_A2_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (with edges)
z_prime_results_JJ_A2 <- grouped_data_JJ_A2 %>%
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
  tab_header(title = "Plate JJ A2 MGM 24h (with edges)")



# edgeless plate
JJ_plate_A2_edgeless <- JJ_plate_A2 %>% 
  filter(!column %in% c("1", "2", "23", "24")) %>% 
  filter(!row %in% c("A", "B", "P", "O")) 
JJ_plate_A2_edgeless_max <- JJ_plate_A2_edgeless %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

# CV for max only (edgeless)
JJ_edgeless_A2_CV_analysis <- JJ_plate_A2_edgeless_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "JJ A2 MGM 24h (without edges)") %>%
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
JJ_A2_Z_data_edgeless <- JJ_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value") %>% 
  # remove outliers
  # filter(!(Treatment == "LPS.dmso" & Threshold == "x2k_nuc" & Value > 1500)) %>%
  group_by(Threshold, Treatment) %>%
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Plot (no edges)
JJ_A2_Z_plot_edgeless <- JJ_A2_Z_data_edgeless %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate JJ A2 MGM 24h (without edges)", y = "CD38 Area/nuclei") 

# Group data by Threshold
grouped_data_JJ_A2_edgeless <- JJ_A2_Z_data_edgeless %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (no edges)
z_prime_results_JJ_A2_edgeless <- grouped_data_JJ_A2_edgeless %>%
  reframe(
    Z_Prime_LPS_TAK3uM = calculate_z_prime(
      Median[Treatment == "LPS.dmso"],
      SD[Treatment == "LPS.dmso"],
      Median[Treatment == "LPS.TAK3mM"],
      SD[Treatment == "LPS.TAK3mM"]
    )
  ) %>% 
  gt() %>% 
  tab_header(title = "Plate JJ A2 MGM 24h (without edges)")


# Baricitinib 1K threshold
JJ_A2_Z_data_Bari <- JJ_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value") 

# Bari Plot
JJ_A2_Z_Bari_plot <- JJ_A2_Z_data_Bari %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Bari 0.3", "LPS.Bari 1", "LPS.Bari 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate JJ A2 (MGM) Baricitinib (24h)", y = "CD38 Area/nuclei")

JJ_A2_Z_data_Bari$Treatment <- relevel(JJ_A2_Z_data_Bari$Treatment, ref = "LPS.dmso")
Bari_lm <- lm(Value ~ factor(Treatment), data = JJ_A2_Z_data_Bari)
plot(Bari_lm)
summary(Bari_lm)
anova(Bari_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Bari 0.3"), c("LPS.dmso", "LPS.Bari 1"), c("LPS.dmso", "LPS.Bari 3"))
JJ_A2_Z_Bari_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Tacro 1K threshold
JJ_A2_Z_data_Tacro <- JJ_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Tacro Plot
JJ_A2_Z_Tacro_plot <- JJ_A2_Z_data_Tacro %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Tac 0.3", "LPS.Tac 1", "LPS.Tac 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate JJ A2 (MGM) Tacrolimus (24h)", y = "CD38 Area/nuclei")

JJ_A2_Z_data_Tacro$Treatment <- relevel(JJ_A2_Z_data_Tacro$Treatment, ref = "LPS.dmso")
Tacro_lm <- lm(Value ~ factor(Treatment), data = JJ_A2_Z_data_Tacro)
plot(Tacro_lm)
summary(Tacro_lm)
anova(Tacro_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Tac 0.3"), c("LPS.dmso", "LPS.Tac 1"), c("LPS.dmso", "LPS.Tac 3"))
JJ_A2_Z_Tacro_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Selinexor 1K threshold
JJ_A2_Z_data_Sel <- JJ_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Sel 3", "Sel 1", "Sel 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Selinexor Plot
JJ_A2_Z_Sel_plot <- JJ_A2_Z_data_Sel %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Sel 0.3", "LPS.Sel 1", "LPS.Sel 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate JJ A2 (MGM) Selinexor (24h)", y = "CD38 Area/nuclei")

JJ_A2_Z_data_Sel$Treatment <- relevel(JJ_A2_Z_data_Sel$Treatment, ref = "LPS.dmso")
Sel_lm <- lm(Value ~ factor(Treatment), data = JJ_A2_Z_data_Sel)
plot(Sel_lm)
summary(Sel_lm)
anova(Sel_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Sel 0.3"), c("LPS.dmso", "LPS.Sel 1"), c("LPS.dmso", "LPS.Sel 3"))
JJ_A2_Z_Sel_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")



# Everolimus 1K threshold
JJ_A2_Z_data_Evero <- JJ_plate_A2_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Ever 3", "Ever 1", "Ever 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,cell_area_1k_rod_h_mglia_dapi_cy5_23) %>% 
  pivot_longer(cols = cell_area_1k_rod_h_mglia_dapi_cy5_23, names_to = "Threshold", values_to = "Value")

# Everolimus Plot
JJ_A2_Z_Evero_plot <- JJ_A2_Z_data_Evero %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Ever 0.3", "LPS.Ever 1", "LPS.Ever 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate JJ A2 (MGM) Everolimus (24h)", y = "CD38 Area/nuclei")

JJ_A2_Z_data_Evero$Treatment <- relevel(JJ_A2_Z_data_Evero$Treatment, ref = "LPS.dmso")
Evero_lm <- lm(Value ~ factor(Treatment), data = JJ_A2_Z_data_Evero)
plot(Evero_lm)
summary(Evero_lm)
anova(Evero_lm)

comparisons_list <- list(c("LPS.dmso", "LPS.Ever 0.3"), c("LPS.dmso", "LPS.Ever 1"), c("LPS.dmso", "LPS.Ever 3"))
JJ_A2_Z_Evero_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")

