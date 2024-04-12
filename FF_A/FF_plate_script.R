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

FF_plate_A1_CM <- FF_plate_A1 %>% 
  filter(media_type == "CM")

FF_plate_A1_MGM <- FF_plate_A1 %>% 
  filter(media_type == "MGM")

# Area analysis
FF_plate_A1_MGM_max <- FF_plate_A1_MGM %>% 
  filter(Compound == "dmso", LPS_status == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
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
JJ_CV_analysis <- JJ_plate_A1_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "JJ A1 CM 24h") %>%
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
# conclusion: use 1k or 2k or 3.5k threshold

# Z data prep (with edges)
FF_A1_MGM_Z_data <- FF_plate_A1_MGM %>% 
  filter(Compound %in% c("dmso", "TAK3uM")) %>%
  # filter(lps_status == c("LPS", "noLPS")) %>% 
  mutate(Treatment = interaction(lps_status, compound)) %>%
  dplyr::select(Treatment,x2k_nuc) %>% 
  pivot_longer(cols = c(x2k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>%
  filter(!Treatment %in% "NoLPS.TAK3uM") 
  summarise(Median = median(Value),
            SD = sd(Value)) 

# JJ_A1_Z_data <- JJ_A1_Z_data %>% 
#   filter(!Treatment %in% "noLPS.dmso")

# Plot (with edges)
FF_A1_MGM_Z_plot <- FF_A1_MGM_Z_data %>% 
  mutate(Treatment = factor(Treatment, levels = c("NoLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot(colour = c("#00bb38", "#f88a82", "#72a5ff")) +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 1) +
  scale_color_manual(values = c(NoLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3uM='#72a5ff')) +
  facet_grid(~ Threshold, labeller = labeller(Threshold = c("x2k_nuc" = "2k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13)) +
  labs(title ="Plate FF A1 MGM 24h (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
  # coord_cartesian(ylim = c(0, 1000))
  
comparisons_list <- list(c("NoLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.TAK3uM"), c("NoLPS.dmso", "LPS.TAK3uM"))
FF_A1_MGM_Z_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")
