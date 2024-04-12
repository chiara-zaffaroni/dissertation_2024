
summary_FF <- EE_7.5K_area_AP %>%
  group_by(lps, compound, Treatment = " ") %>%
  summarise(
    N = n(),
    Mean = mean(area7_5k_nuc),
    Median = median(area7_5k_nuc),
    # Min = min(area7_5k_nuc),
    # Max = max(area7_5k_nuc),
    SD = sd(area7_5k_nuc),
    SEM = sd(area7_5k_nuc) / sqrt(length(area7_5k_nuc)),
    CV = sd(area7_5k_nuc) / mean(area7_5k_nuc) * 100
  ) %>% 
  gt(rowname_col = "Treatment") %>% 
  tab_header(title = "Analysis: Outermost (A+P) Edges") %>%
  fmt_number(decimals = 2) %>% 
  tab_style(
    style = cell_text(
      size = "smaller",
      weight = "bold",
    ),
    locations = cells_body(columns = "CV")
  ) %>% 
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(columns = "CV")
  )