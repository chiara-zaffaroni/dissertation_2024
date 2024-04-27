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
install.packages("sjPlot")
library(sjPlot)

# UPLOAD ALL BATCH DATA
FF_batch <- read_excel("final_hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx")
GG_batch <- read_excel("final_hMglia GG A1 HighValue4 Widerange 10 threshold.xlsx")
HH_batch <- read_excel("final_hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx")
II_batch <- read_excel("final_hMglia II A1 MGM BariTacro Focus 24hrs.xlsx")
JJ_batch <- read_excel("final_hMglia JJ A2.xlsx")
KK_batch <- read_excel("final_hMglia KK A1.xlsx")
LL_batch <- read_excel("final_hMglia LL A1 BariTacro Round4 15.3.24.xlsx")

FF_batch_CM <- read_excel("final_hMGlia FF A1 to 3 Analysed Data Area and Nuclei.xlsx", sheet = 2)
GG_batch_CM <- read_excel("final_hMglia GG A1 HighValue4 Widerange 10 threshold.xlsx", sheet = 2)
HH_batch_CM <- read_excel("final_hMglia HH A3 MGM Highvalue4  AA CD38 10 thresholds.xlsx", sheet = 2)
II_batch_CM <- read_excel("final_hMglia II A1 MGM BariTacro Focus 24hrs.xlsx", sheet = 2)
JJ_batch_CM <- read_excel("final_hMglia JJ A2.xlsx", sheet = 2)



all_batches <- bind_rows(FF_batch, GG_batch, HH_batch, II_batch, JJ_batch, KK_batch, LL_batch)

all_batches_CM <- bind_rows(FF_batch_CM, GG_batch_CM, HH_batch_CM, II_batch_CM, JJ_batch_CM)

all_batches_edgeless <- all_batches %>% 
  filter(!Column %in% c("1", "2", "23", "24")) %>% 
  filter(!Row %in% c("A", "B", "P", "O"))  
  
# Z data prep (with edges)
all_batches_Z_data <- all_batches_CM %>% 
  filter(Compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,x2k_nuc, batch) %>% 
  pivot_longer(cols = c(x2k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>% 
  filter(!Treatment %in% "noLPS.TAK3uM") %>% 
  # filter out anomalies
  filter(!Value > 150000)

# Z data prep (ONLY MAX SIGNAL)
all_batches_Z_data_max <- all_batches %>% 
  filter(Compound %in% c("dmso", "TAK3uM")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment,x2k_nuc, batch) %>% 
  pivot_longer(cols = c(x2k_nuc), names_to = "Threshold", values_to = "Value") %>% 
  group_by(Threshold, Treatment) %>% 
  filter(Treatment %in% "LPS.dmso")

# Compute counts (n) for each treatment
sample_counts <- all_batches_Z_data %>%
  group_by(Treatment) %>%
  summarise(n = n())

# Plot (with edges)
all_batches_Z_plot <- all_batches_Z_data %>%
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>%
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot(colour = c("#00bb38", "#f88a82", "#72a5ff")) +
  theme_minimal() +
  geom_jitter(aes(color = Treatment),position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  scale_color_manual(values = c(noLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3uM='#72a5ff')) +
  # facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30" = "3k Threshold", "cell_area_5k_rod_h_mglia_dapi_cy5_33" = "5k Threshold", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37" = "7.5k Threshold", "cell_area10k_rod_h_mglia_dapi_cy5_41" = "10k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(title ="All Batches MGM (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

comparisons_list <- list(c("noLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.TAK3uM"), c("noLPS.dmso", "LPS.TAK3uM"))
all_batches_Z_plot + geom_signif(comparisons = comparisons_list, map_signif_level = FALSE, step_increase = 0.1, colour = "black") + geom_text(data = sample_counts, aes(label = paste("n =", n), x = Treatment, y = -Inf), vjust = 0, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

# # BY BATCH
# take treatment average for each batch:
# group by batch and within each batch group by treatment and take average

all_batches_Z_data_average <- all_batches_Z_data %>% 
  group_by(Treatment, batch) %>% 
  summarise("Average" = mean(Value))

# Define the new batch names
new_batch_names <- c("Batch 1", "Batch 2", "Batch 3", "Batch 4", "Batch 5")

# Change the levels of the 'batch' variable in your dataframe
all_batches_Z_data_average$batch <- factor(all_batches_Z_data_average$batch, levels = unique(all_batches_Z_data_average$batch), labels = new_batch_names)

# Plot
all_batches_Z_plot <- all_batches_Z_data_average %>%
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>%
  ggplot(aes(x=Treatment, y=Average, fill = batch)) +
  geom_col(position = "dodge") + 
  theme_minimal() +
  # geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  # scale_color_manual(values = c(noLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3uM='#72a5ff')) +
  # facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30" = "3k Threshold", "cell_area_5k_rod_h_mglia_dapi_cy5_33" = "5k Threshold", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37" = "7.5k Threshold", "cell_area10k_rod_h_mglia_dapi_cy5_41" = "10k Threshold"))) +
  theme(strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(title ="All Batches (with edges)", y = "Average of CD38 Area/nuclei", fill = "Batch") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

comparisons_list <- list(c("noLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.TAK3uM"), c("noLPS.dmso", "LPS.TAK3uM"))
all_batches_Z_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black")


# Fit linear model with Treatment and batch as predictors
model <- lm(Average ~ factor(Treatment) * batch, data = all_batches_Z_data_average)
summary(model)

# Perform ANOVA to test for interaction effect
anova_result <- anova(model)



GLM <- lm(Value ~ factor(Treatment), data = all_batches_Z_data)
plot(GLM) # ASSUMPTIONS MET
summary(GLM)
anova(GLM) #F-statistic: 548.2 on 6 and 849 DF,  p-value: < 2.2e-16 --> STRONG EVIDENCE FOR AN EFFECT OF BATCH ON VALUE AKA BATCH-TO-BATCH VARIATION
MEANS <- emmeans(GLM, "Treatment")
PAIRS <- pairs(MEANS)
CI <- confint(PAIRS)

tab_model(GLM, show.se = TRUE, show.std = TRUE, show.stat = TRUE)

tab_model(MEANS)

# Define a vector of new batch names corresponding to the existing levels
new_batch_names <- c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6", "Batch7")

# Change the levels of the 'batch' variable in the data frame
all_batches_Z_data_max$batch <- factor(all_batches_Z_data_max$batch, levels = unique(all_batches_Z_data_max$batch), labels = new_batch_names)
comparisons_list <- combn(new_batch_names, 2, simplify = FALSE)


# Create a function to compute significance and filter non-significant comparisons
compute_significance <- function(comparisons_list, data) {
  # List to store significant comparisons
  significant_comparisons <- list()
  
  # Loop through each comparison
  for (pair in comparisons_list) {
    # Extract batch names
    batch1 <- pair[1]
    batch2 <- pair[2]
    
    # Subset data for current batch pair
    batch_data <- data[data$batch %in% c(batch1, batch2), ]
    
    # Perform statistical test (e.g., t-test) and check significance
    p_value <- t.test(Value ~ batch, data = batch_data)$p.value
    
    # p_value <- PAIRS$p.value
    
    # Check significance level (e.g., p-value threshold for NS)
    if (p_value > 0.05) {  # Adjust threshold as needed
      significant_comparisons[[paste(batch1, "-", batch2)]] <- pair
    }
  }
  
  return(significant_comparisons)
}

# Filter comparisons to show only non-significant (NS) pairs
significant_comparisons <- compute_significance(comparisons_list, all_batches_Z_data)

# ______________
# Compute batch means using emmeans
batch_means <-  emmeans::emmeans(GLM, ~ batch)  # Replace 'model' with your actual model formula

# Perform pairwise comparisons with adjusted p-values
pairwise_comparisons <- pairs(batch_means, adjust = "tukey")  # Adjust method as needed

# Filter non-significant (NS) comparisons based on desired significance level
significant_comparisons <- pairwise_comparisons$prob > 0.05  # Adjust threshold as needed

# Convert significant comparisons to a list format (similar to your existing approach)
significant_comparisons_list <- lapply(which(significant_comparisons), function(i) {
  pair <- pairwise_comparisons$contrasts[i, c("batch1", "batch2")]
  paste(pair, collapse = "-")
})


all_batches_Z_plot <- all_batches_Z_data_max %>% 
  # mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>%
  ggplot(aes(x=batch, y=Value)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5, colour = "darkblue") +
  # scale_color_manual(values = c(noLPS.dmso='#00bb38',LPS.dmso='#f88a82',LPS.TAK3uM='#72a5ff')) +
  # facet_grid(~ Threshold, labeller = labeller(Threshold = c("cell_area_1k_rod_h_mglia_dapi_cy5_23" = "1k Threshold", "cell_area_2k_rod_h_mglia_dapi_cy5_26" = "2k Threshold", "cell_area_3_5k_rod_h_mglia_dapi_cy5_30" = "3k Threshold", "cell_area_5k_rod_h_mglia_dapi_cy5_33" = "5k Threshold", "cell_area_7_5k_rod_h_mglia_dapi_cy5_37" = "7.5k Threshold", "cell_area10k_rod_h_mglia_dapi_cy5_41" = "10k Threshold"))) +
  theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold"), axis.title=element_text(size=13), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(title ="All Batches (with edges)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20)) 

# Compute counts (n) for each batch
batch_counts <- all_batches_Z_data_max %>%
  group_by(batch) %>%
  summarise(n = n())


all_batches_Z_plot + geom_signif(comparisons = significant_comparisons_list, annotation = "NS", map_signif_level = FALSE, step_increase = 0.1, colour = "black") + geom_text(data = batch_counts, aes(label = paste("n =", n), x = batch, y = -Inf), vjust = -0.1, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

# # violin plot
# all_batches_Z_plot <- all_batches_Z_data %>% 
#   mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.TAK3uM"))) %>% 
#   ggplot(aes(x = Treatment, y = Value, fill = Treatment)) +  # Use fill instead of colour
#   geom_violin(colour = "black", alpha = 0.5) +  # Set outline color and transparency
#   geom_jitter(aes(color = Treatment), position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
#   scale_fill_manual(values = c("noLPS.dmso" = '#00bb38', "LPS.dmso" = '#f88a82', "LPS.TAK3uM" = '#72a5ff')) +
#   scale_color_manual(values = c("noLPS.dmso" = '#00bb38', "LPS.dmso" = '#f88a82', "LPS.TAK3uM" = '#72a5ff')) +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 12, face = "bold"),
#         axis.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14)) +
#   labs(title = "All Batches (with edges)", y = "CD38 Area/nuclei", fill = "Treatment Type", color = "Treatment Type") +
#   theme(axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20))

# DRUGS
# Baricitinib 2K threshold
all_batches_Z_data_Bari <- all_batches_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3uM", "Bari 3", "Bari 1", "Bari 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, x2k_nuc, batch) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value") %>% 
  filter(!Treatment %in% c("noLPS.TAK3uM", "noLPS.Bari 0.3", "noLPS.Bari 1", "noLPS.Bari 3"))

# Compute counts (n) for each treatment
sample_counts <- all_batches_Z_data_Bari %>%
  group_by(Treatment) %>%
  summarise(n = n())

# Bari Plot
all_batches_Z_Bari_plot <- all_batches_Z_data_Bari %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.Bari 0.3", "LPS.Bari 1", "LPS.Bari 3", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Baricitinib (all batches)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) 

comparisons_list <- list(c("noLPS.dmso", "LPS.dmso"), c("LPS.dmso", "LPS.Bari 0.3"), c("LPS.dmso", "LPS.Bari 1"), c("LPS.dmso", "LPS.Bari 3"), c("LPS.Bari 0.3", "LPS.TAK3uM"), c("LPS.Bari 1", "LPS.TAK3uM"), c("LPS.Bari 3", "LPS.TAK3uM"))
all_batches_Z_Bari_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black") + geom_text(data = sample_counts, aes(label = paste("n =", n), x = Treatment, y = -Inf), vjust = 0, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

all_batches_Z_data_Bari$Treatment <- relevel(all_batches_Z_data_Bari$Treatment, ref = "LPS.dmso")
GLM <- lm(Value ~ factor(Treatment), data = all_batches_Z_data_Bari)

plot(GLM) # ASSUMPTIONS MET
summary(GLM)
anova(GLM) #F-statistic: 548.2 on 6 and 849 DF,  p-value: < 2.2e-16 --> STRONG EVIDENCE FOR AN EFFECT OF BATCH ON VALUE AKA BATCH-TO-BATCH VARIATION
MEANS <- emmeans(GLM, "Treatment")
PAIRS <- pairs(MEANS)
CI <- confint(PAIRS)


# Tacro 2K threshold
all_batches_Z_data_Tacro <- all_batches_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3uM", "Tac 3", "Tac 1", "Tac 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, x2k_nuc, batch) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value") %>% 
  filter(!Treatment %in% c("noLPS.TAK3uM", "noLPS.Tac 0.3", "noLPS.Tac 1", "noLPS.Tac 3"))

# Compute counts (n) for each treatment
sample_counts <- all_batches_Z_data_Tacro %>%
  group_by(Treatment) %>%
  summarise(n = n())

# Tacro Plot
all_batches_Z_Tacro_plot <- all_batches_Z_data_Tacro %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.Tac 0.3", "LPS.Tac 1", "LPS.Tac 3", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Tacrolimus (all batches)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) 

comparisons_list <- list(c("LPS.dmso", "LPS.Tac 0.3"), c("LPS.dmso", "LPS.Tac 1"), c("LPS.dmso", "LPS.Tac 3"), c("LPS.Tac 0.3", "LPS.TAK3uM"), c("LPS.Tac 1", "LPS.TAK3uM"), c("LPS.Tac 3", "LPS.TAK3uM"))
all_batches_Z_Tacro_plot + geom_signif(comparisons = comparisons_list, map_signif_level = FALSE, step_increase = 0.1, colour = "black") + geom_text(data = sample_counts, aes(label = paste("n =", n), x = Treatment, y = -Inf), vjust = 0, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

all_batches_Z_data_Tacro$Treatment <- relevel(all_batches_Z_data_Tacro$Treatment, ref = "LPS.dmso")
GLM <- lm(Value ~ factor(Treatment), data = all_batches_Z_data_Tacro)

plot(GLM) # ASSUMPTIONS MET
summary(GLM)
# Selinexor 2K threshold
all_batches_Z_data_Sel <- all_batches_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3uM", "Sel 3", "Sel 1", "Sel 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, x2k_nuc, batch) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value") %>% 
  filter(!Treatment %in% c("noLPS.TAK3uM", "noLPS.Sel 0.3", "noLPS.Sel 1", "noLPS.Sel 3"))

# Compute counts (n) for each treatment
sample_counts <- all_batches_Z_data_Sel %>%
  group_by(Treatment) %>%
  summarise(n = n())

# Selinexor Plot
all_batches_Z_Sel_plot <- all_batches_Z_data_Sel %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.Sel 0.3", "LPS.Sel 1", "LPS.Sel 3", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Selinexor (all batches)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) 

comparisons_list <- list(c("LPS.dmso", "LPS.Sel 0.3"), c("LPS.dmso", "LPS.Sel 1"), c("LPS.dmso", "LPS.Sel 3"), c("LPS.Sel 0.3", "LPS.TAK3uM"), c("LPS.Sel 1", "LPS.TAK3uM"), c("LPS.Sel 3", "LPS.TAK3uM"))
all_batches_Z_Sel_plot + geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, step_increase = 0.1, colour = "black") + geom_text(data = sample_counts, aes(label = paste("n =", n), x = Treatment, y = -Inf), vjust = 0, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

all_batches_Z_data_Sel$Treatment <- relevel(all_batches_Z_data_Sel$Treatment, ref = "LPS.dmso")
GLM <- lm(Value ~ factor(Treatment), data = all_batches_Z_data_Sel)
summary(GLM)

# Everolimus 2K threshold
all_batches_Z_data_Evero <- all_batches_edgeless %>% 
  filter(Compound %in% c("dmso", "TAK3uM", "Ever 3", "Ever 1", "Ever 0.3")) %>% 
  mutate(Treatment = interaction(LPS_status, Compound)) %>% 
  dplyr::select(Treatment, x2k_nuc, batch) %>% 
  pivot_longer(cols = x2k_nuc, names_to = "Threshold", values_to = "Value") %>% 
  filter(!Treatment %in% c("noLPS.TAK3uM", "noLPS.Ever 0.3", "noLPS.Ever 1", "noLPS.Ever 3"))

# Compute counts (n) for each treatment
sample_counts <- all_batches_Z_data_Evero %>%
  group_by(Treatment) %>%
  summarise(n = n())

# Everolimus Plot
all_batches_Z_Evero_plot <- all_batches_Z_data_Evero %>% 
  mutate(Treatment = factor(Treatment, levels = c("noLPS.dmso", "LPS.dmso", "LPS.Ever 0.3", "LPS.Ever 1", "LPS.Ever 3", "LPS.TAK3uM"))) %>% 
  ggplot(aes(x=Treatment, y=Value, colour = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  facet_wrap(~ Threshold) + 
  theme(legend.position = "none")+
  labs(title ="Everolimus (all batches)", y = "CD38 Area/nuclei") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) 

comparisons_list <- list(c("LPS.dmso", "LPS.Ever 0.3"), c("LPS.dmso", "LPS.Ever 1"), c("LPS.dmso", "LPS.Ever 3"))
all_batches_Z_Evero_plot + geom_signif(comparisons = comparisons_list, map_signif_level = FALSE, step_increase = 0.1, colour = "black") + geom_text(data = sample_counts, aes(label = paste("n =", n), x = Treatment, y = -Inf), vjust = 0, size = 4, hjust = 0.5, angle = 0, inherit.aes = FALSE)

all_batches_Z_data_Evero$Treatment <- relevel(all_batches_Z_data_Evero$Treatment, ref = "LPS.dmso")
GLM <- lm(Value ~ factor(Treatment), data = all_batches_Z_data_Evero)
summary(GLM)
