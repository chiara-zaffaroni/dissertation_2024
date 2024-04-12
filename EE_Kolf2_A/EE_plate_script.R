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
  clean_names() %>% 
  filter(media_type %in% "MGM") %>% 
  select(media_type, lps, compound, area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc)

EE_plate_max <- EE_plate %>% 
  filter(compound == "dmso", lps == "LPS") %>% 
  pivot_longer(cols = contains("area"), names_to = "Variable", values_to = "Value") %>% 
  select(Variable, Value) %>% 
  mutate(Value = as.numeric(Value))

EE_plate_max_plot <- EE_plate_max %>% 
  # edit naming of x axis tick labels and arrange in numerical order
  mutate(x_label = gsub("^cell_area10k.*", "cell_area10k", Variable)) %>%
  mutate(x_label = gsub("^(cell_area_[0-9]+(_[0-9]+)?k).*", "\\1", x_label)) %>%  
  ggplot(aes(x = fct_inorder(x_label), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Plate EE (MGM)",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")

# CV for max only (with edges)
EE_CV_analysis <- EE_plate_max %>% 
  group_by(Variable) %>% 
  # select(Variable, Value) %>% 
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "EE MGM 24h") %>%
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




# remove edges (1+2, 23+24, A+B, O+P)
EE_plate_edgefree <- EE_plate %>% 
  filter(!column %in% c("1", "2", "23", "24")) %>% 
  filter(!row %in% c("A", "B", "P", "O"))

# 2. MIN/MAX ALL THRESHOLDS EDGELESS
# separate by media type and area
# CM + area
EE_plate_edgefree_CM_area <- EE_plate_edgefree %>% 
  filter(media_type %in% "CM") %>% 
  select(media_type, lps, compound, area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc)

# CM + LPS/dmso area
EE_plate_edgefree_CM_area_LPS_dmso <- EE_plate_edgefree_CM_area %>% 
  filter(lps %in% "LPS" & compound %in% "dmso")

# CM + noLPS/dmso area
EE_plate_edgefree_CM_area_noLPS_dmso <- EE_plate_edgefree_CM_area %>% 
  filter(lps %in% "NoLPS" & compound %in% "dmso")

# plotting: CM + LPS/dmso area
EE_plate_edgefree_CM_area_LPS_dmso_long <- EE_plate_edgefree_CM_area_LPS_dmso %>% 
  pivot_longer(cols = c(area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc), names_to = "Variable", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value))

area_CM_LPS_dmso_plot <- EE_plate_edgefree_CM_area_LPS_dmso_long %>% 
  ggplot(aes(x = factor(Variable, levels = c("area5k_nuc", "area7_5k_nuc", "area10k_nuc", "area12_5k_nuc", "area15k_nuc")), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, colour = "darkblue") +
  labs(title = "Maximun Signal (LPS/dmso) in CM",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")

# plotting: CM + noLPS/dmso area
EE_plate_edgefree_CM_area_noLPS_dmso_long <- EE_plate_edgefree_CM_area_noLPS_dmso %>% 
  pivot_longer(cols = c(area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc), names_to = "Variable", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value))

area_CM_noLPS_dmso_plot <- EE_plate_edgefree_CM_area_noLPS_dmso_long %>% 
  ggplot(aes(x = factor(Variable, levels = c("area5k_nuc", "area7_5k_nuc", "area10k_nuc", "area12_5k_nuc", "area15k_nuc")), y = Value)) +
  geom_boxplot(aes(colour = "darkblue")) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "Minimum Signal (noLPS/dmso) in CM",
       x = "Threshold",
       y = "CD38 Area/nuclei") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = "darkblue")

# Find the maximum value across both plots
max_value <- max(c(max(EE_plate_edgefree_CM_area_LPS_dmso_long$Value), max(EE_plate_edgefree_CM_area_noLPS_dmso_long$Value)))

# Set the same y-axis limits for both plots
area_CM_LPS_dmso_plot <- area_CM_LPS_dmso_plot + coord_cartesian(ylim = c(0, max_value))
area_CM_noLPS_dmso_plot <- area_CM_noLPS_dmso_plot + coord_cartesian(ylim = c(0, max_value))

grid.arrange(area_CM_LPS_dmso_plot, area_CM_noLPS_dmso_plot, ncol=2)



# CV: CM + LPS/dmso area 
sapply()

# LPS/TAK3 area
# NoLPS/dmso area




# 3. MEDIA EFFECT
EE_plate_area_LPS_dmso <- EE_plate_edgefree %>%
  # filter for lps/compound combo first
  filter(lps == "LPS" & compound == "dmso") %>% 
  # group by media type
  group_by(media_type) %>% 
  select(media_type, lps, compound, area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc)

EE_plate_area_LPS_dmso_long <- EE_plate_area_LPS_dmso %>% 
  pivot_longer(cols = c(area5k_nuc, area7_5k_nuc, area10k_nuc, area12_5k_nuc, area15k_nuc), names_to = "Variable", values_to = "Value") %>% 
  mutate(Value = as.numeric(Value))

EE_plate_area_LPS_dmso_plot <- EE_plate_area_LPS_dmso_long %>% 
  ggplot(aes(x = factor(Variable, levels = c("area5k_nuc", "area7_5k_nuc", "area10k_nuc", "area12_5k_nuc", "area15k_nuc")), y = Value, colour = media_type)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3) +
  labs(title = "Area all LPS/dmso (Max)",
       x = "Threshold",
       y = "CD38 Area/nuclei",
       colour = "Media Type") +
  theme_minimal() +
  scale_color_manual(values = c("darkblue", "darkgreen"))

