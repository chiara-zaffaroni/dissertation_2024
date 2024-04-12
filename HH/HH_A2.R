HH_plate_A2_original <- read_excel("hMglia HH A2 CM MinMax LPSdmso vs LPS TAK242 3uM.xlsx", sheet = 1, skip = 8) %>% 
  clean_names() %>% 
  # separate row and column
  separate(well_name, into = c("Row", "Column"), sep = -2) %>% 
  mutate(Column = as.numeric(Column)) %>% 
  arrange(Column) %>% 
  # add media type column
  mutate(Media_Type = "CM") %>% 
  select(Media_Type, Row, Column, nuclei_count_total_rod_h_mglia_dapi_cy5_14, everything()) %>% 
  rename("Nuclei_Count" ="nuclei_count_total_rod_h_mglia_dapi_cy5_14") %>%
  # select_if(~ !any(is.na(.))) %>%
  select(-c(plate_id, measurement_set_id, cell_count_rod_h_mglia_dapi_cy5_4, cell_object_id_rod_h_mglia_dapi_cy5_5,
            cell_10k_objects_rod_h_mglia_dapi_cy5_13, cell_10k_objects_rod_h_mglia_dapi_cy5_72,cell_12_5k_objects_rod_h_mglia_dapi_cy5_106,
            cell_12_5k_objects_rod_h_mglia_dapi_cy5_47, cell_15k_objects_rod_h_mglia_dapi_cy5_66,
            cell_15k_objects_rod_h_mglia_dapi_cy5_7, cell_17_5k_objects_rod_h_mglia_dapi_cy5_113,
            cell_17_5k_objects_rod_h_mglia_dapi_cy5_54, cell_1k_objects_rod_h_mglia_dapi_cy5_11, 
            cell_1k_objects_rod_h_mglia_dapi_cy5_70, cell_20k_objects_rod_h_mglia_dapi_cy5_117,
            cell_20k_objects_rod_h_mglia_dapi_cy5_58, cell_25k_objects_rod_h_mglia_dapi_cy5_121,
            cell_25k_objects_rod_h_mglia_dapi_cy5_62, cell_2k_objects_rod_h_mglia_dapi_cy5_29,
            cell_2k_objects_rod_h_mglia_dapi_cy5_88, cell_3_5k_objects_rod_h_mglia_dapi_cy5_12,
            cell_3_5k_objects_rod_h_mglia_dapi_cy5_71, cell_5k_objects_rod_h_mglia_dapi_cy5_36,
            cell_5k_objects_rod_h_mglia_dapi_cy5_95, cell_7_5k_objects_rod_h_mglia_dapi_cy5_40,
            cell_7_5k_objects_rod_h_mglia_dapi_cy5_99, cell_count_rod_h_mglia_dapi_cy5_63,
            cell_nuclei_count_rod_h_mglia_dapi_cy5_22, cell_object_id_rod_h_mglia_dapi_cy5_64,
            x15k_objects_total_rod_h_mglia_dapi_cy5_65, x3_5k_objects_total_rod_h_mglia_dapi_cy5_68)) %>%
  select(1:48) %>% 
  mutate_at(vars(5:48), ~ . / Nuclei_Count)


HH_LPS_map_A2 <- read_excel("hMglia HH A2 CM MinMax LPSdmso vs LPS TAK242 3uM.xlsx", sheet = 3) %>% 
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Value", everything())

HH_platemap_A2 <- read_excel("hMglia HH A2 CM MinMax LPSdmso vs LPS TAK242 3uM.xlsx", sheet = 2, skip = 4) %>% 
  rename(Row = "...1") %>%
  mutate(row_names = .[[1]]) %>%
  select(-1) %>%
  column_to_rownames(var = "row_names") %>% 
  gather(key = "Column", value = "Compound", everything()) %>% 
  bind_cols(HH_LPS_map_A2$Value) %>% 
  rename(LPS_status = "...3")



HH_plate_A2 <- HH_plate_A2_original %>% 
  bind_cols(HH_platemap_A2[c("Compound", "LPS_status")]) %>% 
  # rename(Compound = "...49") %>%
  select(Media_Type, Row, Column, Compound, LPS_status, Nuclei_Count, everything())

# export into an excel doc
write_xlsx(HH_plate_A2, "final_hMglia HH A2 CM MinMax LPSdmso vs LPS TAK242 3uM.xlsx")
