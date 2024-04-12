GG_plate_A2_original <- read_excel("hMglia GG A2 MinMax LPS vs No LPS.xlsx", sheet = 1, skip = 8) %>% 
  clean_names() %>% 
  # separate row and column
  separate(well_name, into = c("Row", "Column"), sep = -2) %>% 
  mutate(Column = as.numeric(Column)) %>% 
  arrange(Column) %>% 
  # add media type column
  mutate(Media_Type = case_when(
    row_number() <= 192 ~ "CM",
    row_number() > 192 ~ "MGM"
  )) %>% 
  select(Media_Type, Row, Column, cell_nuclei_count_rod_h_mglia_dapi_cy5, everything()) %>% 
  rename("Nuclei_Count" ="cell_nuclei_count_rod_h_mglia_dapi_cy5") %>%
  select_if(~ !any(is.na(.))) %>%
  select(-c(plate_id, measurement_set_id, cell_object_id_rod_h_mglia_dapi_cy5)) %>%
  select(1:48) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)


GG_LPS_map_A2 <- read_excel("hMglia GG A2 MinMax LPS vs No LPS.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

GG_platemap_A2 <- read_excel("hMglia GG A2 MinMax LPS vs No LPS.xlsx", sheet = 2, skip = 5)

GG_platemap_A2 <- GG_platemap_A2[1:(nrow(GG_platemap_A2) - 4), ] %>% 
  # rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(GG_LPS_map_A2$Value) %>% 
  rename(LPS_status = "...3")



GG_plate_A2 <- GG_plate_A2_original %>% 
  bind_cols(GG_platemap_A2[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(GG_plate_A2, "final_hMglia GG A2 MinMax LPS vs No LPS.xlsx")
