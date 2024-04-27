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
library(emmeans)
# read file
EE_plate <-  read_excel("hMglia EE OM Good Plate Kolf2 A.xlsx", sheet=1, skip=5) %>% 
  clean_names()

# Get the first row of the data as a vector
new_headers <- as.character(EE_plate[1, ])

# Replace the existing column names with the new headers
colnames(EE_plate) <- new_headers

# Remove the first row from the data frame
EE_plate <- EE_plate[-1, ] %>% 
  clean_names()



FF_plate_A1_original <- read_excel("hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx", sheet = 1, skip = 10) %>% 
  clean_names() %>% 
  rename("Nuclei_Count" ="nuclei_count_total_rod_h_mglia_dapi_cy5") %>% 
  select(media_type, row, column, lps_status, compound, Nuclei_Count, x2k_nuc, x3_5k_nuc, x5k_nuc, x7_5k_nuc, x10k_nuc, x12_5k_nuc, x15k_nuc, x17_5k_nuc, x20k_nuc, x25k_nuc)

FF_LPS_map_A1 <- read_excel("hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx", sheet = 5) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

FF_platemap_A1 <- read_excel("hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx", sheet = 4, skip = 10) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(FF_LPS_map_A1$Value) %>% 
  rename(LPS_status = "...3")

FF_plate_A1 <- FF_plate_A1_original %>% 
  bind_cols(FF_platemap_A1[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(media_type, row, column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(FF_plate_A1, "final_hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx")     

FF_plate_A1 <- read_excel("final_hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx")

FF_plate_A1_CM <- FF_plate_A1 %>% 
  filter(media_type == "CM")

FF_plate_A1_MGM <- FF_plate_A1 %>% 
  filter(media_type == "MGM") %>% 
  select(media_type, everything())

FF_plate_A1_MGM_edgeless <- FF_plate_A1 %>% 
  filter(media_type == "MGM") %>% 
  filter(!column %in% c("1", "2", "23", "24")) %>% 
  filter(!row %in% c("A", "B", "P", "O")) %>% 
  select(media_type, everything())

# Area analysis
FF_plate_A1_MGM_max <- FF_plate_A1_MGM %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("k_nuc"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

FF_plate_A1_MGM_edgeless_max <- FF_plate_A1_MGM_edgeless %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("k_nuc"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value)

FF_plate_A1_area_plot <- FF_plate_A1_MGM_max %>% 
  # edit naming of x axis tick labels and arrange in numerical order
  mutate(x_label = gsub("^cell_area10k.*", "cell_area10k", Variable)) %>%
  mutate(x_label = gsub("^(cell_area_[0-9]+(_[0-9]+)?k).*", "\\1", x_label)) %>%  
  ggplot(aes(x = fct_inorder(x_label), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Plate FF A1 (CM 24h)",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")

# CV for max only (with edges)
FF_MGM_CV_analysis <- FF_plate_A1_MGM_edgeless_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "FF A1 MGM (without edges)") %>%
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

# Z data prep (with edges)
FF_A1_MGM_Z_data <- FF_plate_A1_MGM %>% 
  filter(Compound %in% c("dmso", "TAK3uM")) %>%
  # filter(lps_status == c("LPS", "noLPS")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>%
  dplyr::select(Treatment,x2k_nuc) %>% 
  pivot_longer(cols = c(x2k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  filter(!Treatment %in% "noLPS.TAK3uM") %>%  
  summarise(Median = median(Value),
            SD = sd(Value)) 

# Compute counts (n) for each treatment
sample_counts <- FF_A1_MGM_Z_data %>%
  group_by(Treatment) %>%
  summarise(n = n())

# Plot (with edges)
FF_A1_MGM_Z_plot <- FF_A1_MGM_Z_data %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot(colour = c("#00bb38", "#f88a82", "#72a5ff")) +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  scale_color_manual(values = c(noLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3uM='#72a5ff')) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("x2k_nuc" = "2k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(title ="Plate FF A1 MGM 24h (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 3500))
  
comparisons_list <- list(c("noLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.TAK3uM"), c("noLPS.dmso", "LPS.TAK3uM"))
FF_A1_MGM_Z_plot + geom_signif(comparisons = comparisons_list, map_signif_level = FALSE, step_increase = 0.1, colour = "black") + geom_text(data = sample_counts, aes(label = paste("n =", n), x = Treatment, y = -Inf), vjust = 0, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

GLM <- lm(Value ~ factor(Treatment), data = FF_A1_MGM_Z_data)
plot(GLM) # ASSUMPTIONS MET
summary(GLM)
anova(GLM) #F-statistic: 548.2 on 6 and 849 DF,  p-value: < 2.2e-16 --> STRONG EVIDENCE FOR AN EFFECT OF BATCH ON VALUE AKA BATCH-TO-BATCH VARIATION
MEANS <- emmeans(GLM, "Treatment")
PAIRS <- pairs(MEANS)
CI <- confint(PAIRS)


# Group data by Threshold
grouped_data_FF_A1 <- FF_A1_MGM_Z_data %>%
  group_by(Threshold)

# Define function to calculate Z' prime
calculate_z_prime <- function(positive_median, positive_sd, negative_median, negative_sd) {
  Z_Prime <- 1 - (3 * (negative_sd + positive_sd)) / abs(negative_median - positive_median)
  return(Z_Prime)
}

# Calculate Z' prime for both combinations (no edges)
z_prime_results_FF_A1 <- grouped_data_FF_A1 %>%
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
  tab_header(title = "Plate FF A1 MGM (without edges)")



# Baricitinib 1K threshold
FF_A1_Z_data_Bari <- FF_plate_A1_MGM_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,x2k_nuc) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value") 

# Bari Plot
FF_A1_Z_Bari_plot <- FF_A1_Z_data_Bari %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Bari 0.3", "LPS.Bari 1", "LPS.Bari 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) +
  theme(legend.position = "none")+
  labs(title ="Plate FF A1 Baricitinib", y = "CD38 Area/nuclei")


comparisons_list <- list(c("LPS.dmso", "LPS.Bari 0.3"), c("LPS.dmso", "LPS.Bari 1"), c("LPS.dmso", "LPS.Bari 3"))
FF_A1_Z_Bari_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# # calculating the mean values for each type of Treatment
# Bari_emmeans <- emmeans(Bari_lm, "Treatment")
# Bari_pairs <- pairs(Bari_emmeans)
# confint(Bari_pairs)

# Tacro 1K threshold
FF_A1_Z_data_Tacro <- FF_plate_A1_MGM_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,x2k_nuc) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value")

# Tacro Plot
FF_A1_Z_Tacro_plot <- FF_A1_Z_data_Tacro %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Tac 0.3", "LPS.Tac 1", "LPS.Tac 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate FF A1 Tacrolimus", y = "CD38 Area/nuclei")

comparisons_list <- list(c("LPS.dmso", "LPS.Tac 0.3"), c("LPS.dmso", "LPS.Tac 1"), c("LPS.dmso", "LPS.Tac 3"))
FF_A1_Z_Tacro_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Selinexor 1K threshold
FF_A1_Z_data_Sel <- FF_plate_A1_MGM_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Sel 3", "Sel 1", "Sel 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, x2k_nuc) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value")

# Selinexor Plot
FF_A1_Z_Sel_plot <- FF_A1_Z_data_Sel %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Sel 0.3", "LPS.Sel 1", "LPS.Sel 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate FF A1 Selinexor", y = "CD38 Area/nuclei")
comparisons_list <- list(c("LPS.dmso", "LPS.Sel 0.3"), c("LPS.dmso", "LPS.Sel 1"), c("LPS.dmso", "LPS.Sel 3"))
FF_A1_Z_Sel_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Everolimus 1K threshold
FF_A1_Z_data_Evero <- FF_plate_A1_MGM_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3mM", "Ever 3", "Ever 1", "Ever 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,x2k_nuc) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value")

# Everolimus Plot
FF_A1_Z_Evero_plot <- FF_A1_Z_data_Evero %>% 
  mutate(Treatment = factor(Treatment, levels = c("LPS.dmso", "LPS.Ever 0.3", "LPS.Ever 1", "LPS.Ever 3", "LPS.TAK3mM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Plate FF A1 Everolimus", y = "CD38 Area/nuclei")
comparisons_list <- list(c("LPS.dmso", "LPS.Ever 0.3"), c("LPS.dmso", "LPS.Ever 1"), c("LPS.dmso", "LPS.Ever 3"))
FF_A1_Z_Evero_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")