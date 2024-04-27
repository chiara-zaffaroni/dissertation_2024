data<-read.table("Rescue.csv", header=TRUE, sep=",") %>%
  group_by(Cell.line, Treatment)

# data_long <- data %>%
#   pivot_wider(names_from = Treatment, values_from = No..of.individual.cells) %>%
#   pivot_longer(cols = c(DMSO, FSK), names_to = "Treatment", values_to = "No..of.individual.cells")


datagraph<-data %>%
  # mutate(Cell.line = factor(Cell.line, levels = c("Wild type", "FAM134B1C3+V5", 	
  #                                                 "FAM134B1C3+FAMB", 	
  #                                                 "FAM134B1C3+FAMC", "FAM134B3C1+V5", "FAM134B3C1+FAMB", 	
  #                                                 "FAM134B3C1+FAMC")))+
  ggplot(aes(x=Cell.line, y=No..of.individual.cells, fill = Treatment)+
  geom_bar()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title="Rescue of phenotype", x="Cell line", y="%Number of individual cells")

  