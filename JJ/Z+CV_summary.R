library(reshape2)

CV_data <- read_excel("Summary CV and Z prime for CM vs MGM.xlsx", sheet = 2)

Z_data <- read_excel("Summary CV and Z prime for CM vs MGM.xlsx", sheet = 3)


CV_data_melted <- melt(CV_data)

CV_plot <- CV_data_melted %>% 
  ggplot(aes(x = variable, y = value, colour = variable)) +
  geom_point() +
  theme_minimal() +
  geom_point(position = position_jitter(width = 0.1, height = 0.1), size = 4) + 
  labs(x = "Media Type", y = "Coefficient of Variation (%)") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = c("CM" = "#ff7700", "MGM" = "blue"))



Z_data_melted <- melt(Z_data)

Z_plot <- Z_data_melted %>% 
  ggplot(aes(x = variable, y = value, colour = variable)) +
  geom_point() +
  theme_minimal() +
  geom_point(size = 4) +
  labs(x = "Media Type", y = "Z prime") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(-3, 1)) +
  scale_color_manual(values = c("CM" = "#ff7700", "MGM" = "blue"))
