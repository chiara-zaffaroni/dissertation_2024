library(gtsummary)
library(gt)

# CM
GG_CV_analysis_CM <- GG_A1_max_area_long %>% 
  filter(Media_Type == "CM") %>% 
  group_by(Threshold_using_area) %>% 
  select(Threshold_using_area, Value) %>% 
  # mutate(Mean = mean(Value))
  summarise(
    # N = n(),
    # Mean = mean(Value),
    # Median = median(Value),
    # Min = min(area7_5k_nuc),
    # Max = max(area7_5k_nuc),
    # SD = sd(Value),
    # SEM = sd(Value) / sqrt(length(Value)),
    CV = sd(Value) / mean(Value) * 100
  ) %>% 
  gt() %>% 
  tab_header(title = "CM media") %>%
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
    style = list(cell_text(color = 'lightgreen'))
  )
# conclusion: select 1K, 2K

# MGM
GG_CV_analysis_MGM <- GG_A1_max_area_long %>% 
  filter(Media_Type == "MGM") %>% 
  group_by(Threshold_using_area) %>% 
  select(Threshold_using_area, Value) %>% 
  # mutate(Mean = mean(Value))
  summarise(CV = sd(Value) / mean(Value) * 100) %>% 
  gt() %>% 
  tab_header(title = "MGM media") %>%
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
    style = list(cell_text(color = 'lightgreen'))
  )
# conclusion: select 1K, 2K, 3.5K
# then looked at heat map and concluded 2K to be cleanest