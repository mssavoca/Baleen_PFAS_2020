--------------------------------------------------------------------------------
  # PFAS tables, summary tables, comparison tables, ratio tables across tissue types 
  # Matt Savoca
  # Started on: 7/18/23
--------------------------------------------------------------------------------



View(PFAS_full_data_comb)


# Define the desired order of compounds
compounds_order <- c("PFOA", "PFOS", "FOSA", "7:3 FTCA", "PFDA", 
                     "PFNA", "PFUdA", "PFtrDA", "PFHxA", "PFHxS", "PFDoA")


# summary tables of PFAS data----


summ_pfas_baleen <- PFAS_full_data_comb %>%
  filter(
    Sample_type == "liver",
    Compound %in% c(PFASofInterest),
    #Sample_seq == 1,
    #Prey == "Zooplankton"
    ) %>%
  group_by(ID_code) %>%
  summarise(
    Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
    Median_Conc = median(Conc_Corr_num, na.rm = TRUE),
    SD_Conc = sd(Conc_Corr_num, na.rm = TRUE)
  ) %>%
  ungroup()

median(summ_pfas_baleen$Total_PFAS)
sd(summ_pfas_baleen$Total_PFAS)


Summ_table_CompoundsListed <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type != "gum-baleen interface") %>% 
  group_by(ID_code, Sample_type, Plate_num, Compound) %>% 
  summarise(
    Med_conc = median(Conc_Corr_num, na.rm = TRUE),
    SD_conc = sd(Conc_Corr_num, na.rm = TRUE),
  ) %>% 
  arrange(match(Compound, compounds_order)) %>% 
  ungroup() %>% 
  mutate(Med_SD = paste(round(Med_conc, 2), " (", round(SD_conc, 2), ")", sep = ""))



Summ_table_CompoundsListed_wide <- pivot_wider(
  data = Summ_table_CompoundsListed,
  id_cols = c(ID_code, Sample_type, Plate_num),
  names_from = Compound,
  values_from = Med_SD,
  names_prefix = ""
) %>% 
  arrange(Sample_type)

write.csv(Summ_table_CompoundsListed_wide, "Summ table all compounds.csv", row.names = FALSE)



Overall_summ_table_baleen <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type == "baleen") %>% 
  group_by(ID_code, Sample_type, Plate_num, Compound) %>% 
  summarise(
    Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
  ) %>% 
  ungroup() %>% 
  group_by(ID_code, Sample_type, Plate_num) %>% 
  summarise(
    Med_conc = median(Total_PFAS, na.rm = TRUE),
    SD_conc = sd(Total_PFAS, na.rm = TRUE),
  ) %>% 
  ungroup() %>% 
  arrange(Sample_type) %>% 
  mutate(Med_SD = paste(round(Med_conc, 2), " (", round(SD_conc, 2), ")", sep = "")) %>% 
  select(ID_code, Sample_type, Plate_num, Med_SD)



Overall_summ_table_WT <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         !Sample_type %in% c("baleen", "gum-baleen interface")) %>% 
  group_by(ID_code, Sample_type, Plate_num) %>% 
  summarise(
    Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
    Med_conc = median(Total_PFAS, na.rm = TRUE),
    SD_conc = sd(Total_PFAS, na.rm = TRUE),
  ) %>% 
  ungroup() %>% 
  arrange(Sample_type) %>% 
  mutate(Med_SD = paste(round(Med_conc, 2), " (",round(SD_conc, 2),")", sep = "")) %>% 
  select(ID_code, Sample_type, Plate_num, Med_SD)
       
Overall_summ_table <- rbind(Overall_summ_table_baleen, Overall_summ_table_WT)


Overall_summ_table_wide <- pivot_wider(
  data = Overall_summ_table,
  id_cols = c(ID_code, Plate_num),
  names_from = Sample_type,
  values_from = Med_SD,
  names_prefix = ""
) %>%
  arrange(ID_code, Plate_num)

write.csv(Overall_summ_table_wide, "Summ table top comps comb.csv", row.names = FALSE)




PFAS_for_study_comp <- PFAS_full_data_comb %>% 
  filter(Compound %in% c("PFOA", "PFOS", "FOSA", "7:3 FTCA", "PFDA", 
                         "PFNA", "PFUdA", "PFtrDA", "PFHxA", "PFHxS", "PFDoA"),
         Sample_type == "liver"
  ) %>% 
  mutate(Compound = factor(Compound, levels = compounds_order)) %>% 
  group_by(Species, Compound) %>% 
  summarise(med_val = median(Conc_Corr_num))



#Ratio tables----



# Tissue-PFAS ratios----

PFAS_by_Tiss_prop <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type != "gum-baleen interface",
         ID_code %in% c("IFAW13-158Mn", "IFAW16-227Mn", "IFAW17-182Eg", "IFAW19-287Mn", "IFAW20-009Mn"),
         ifelse(Sample_type == "baleen", if_else(ID_code == "IFAW13-158Mn", 
                                                 Sample_seq == 2, Sample_seq == 1), TRUE)) %>% 
  group_by(ID_code, Sample_type, Plate_num) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE))
  

ratio_df <- PFAS_by_Tiss_prop %>%
  pivot_wider(
    names_from = c(Sample_type, Plate_num),
    values_from = Total_PFAS,
    values_fill = NA
  ) %>% 
  mutate(
    Ratio_baleen1_liver = baleen_1 / liver_NA,
    Ratio_baleen1_skin = baleen_1 / skin_NA,
    Ratio_baleen1_blubber = baleen_1 / blubber_NA,
    Ratio_baleen1_gum = baleen_1 / gum_NA,
    Ratio_baleen2_liver = baleen_2 / liver_NA,
    Ratio_baleen2_skin = baleen_2 / skin_NA,
    Ratio_baleen2_blubber = baleen_2 / blubber_NA,
    Ratio_baleen2_gum = baleen_2 / gum_NA,
  )



# Compound PFAS ratios, across tissues----
SpecPFASratios_df <- PFAS_full_data_comb %>% 
  filter(Compound %in% c("PFOS", "FOSA", "PFUdA"),
         Sample_type != "gum-baleen interface",
         ID_code %in% c("IFAW13-158Mn", "IFAW16-227Mn", "IFAW17-182Eg", "IFAW19-287Mn", "IFAW20-009Mn"),
         ifelse(Sample_seq == "baleen", if_else(ID_code == "IFAW13-158Mn", 
                                                 Sample_seq == 2, Sample_seq == 1), TRUE)) %>% 
  group_by(ID_code, Sample_type, Plate_num) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE))
         





