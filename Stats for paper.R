--------------------------------------------------------------------------------
  # Stats for PFAS-baleen paper
  # Matt Savoca
  # Started on: 7/28/23
--------------------------------------------------------------------------------

  TL_PFAS_idv <- Baleen_PFAS_data_comb %>%
    filter(Compound %in% PFASofInterest,
           Sample_type == "baleen",
           ifelse(ID_code == "IFAW13-158Mn", Sample_seq == 2, Sample_seq == 1)) %>% 
    group_by(ID_code, Plate_num) %>% 
    summarise(
      Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
      ID_code = first(ID_code),
      SciName = abbr_binom(first(SciName)),
      Troph_level = first(Troph_level),
      Prey = first(Prey),
      Coast = first(Coast), 
      Plate_num = first(Plate_num),
      .groups = 'drop'  # Suppress the grouping message
    ) %>%
    mutate(Prey = ifelse(Prey == "Mixed (ZP/fish)", "Zooplankton and fish",
                         ifelse(Prey == "Zooplankton", "Zooplankton-specialist", Prey))) %>%
    mutate(Prey = factor(Prey, levels = c("Zooplankton-specialist", "Zooplankton and fish")))
  
  
  TL_PFAS_idv_summ <- TL_PFAS_idv %>% 
    group_by(Prey) %>% 
    summarise(
      Med_tot_PFAS = median(Total_PFAS),
      SD_tot_PFAS = sd(Total_PFAS)
    )
  
# Perform Kruskal-Wallis test
PFASbyPrey <- kruskal.test(Total_PFAS ~ Prey, data = TL_PFAS_idv)
  
# Print the result
print(PFASbyPrey)





FOSA_baleen_idv <- Baleen_PFAS_data_comb %>%
  filter(Compound %in% PFASofInterest,
         Sample_type == "baleen",
         Compound == "FOSA",
         ifelse(ID_code == "IFAW13-158Mn", Sample_seq == 2, Sample_seq == 1)) %>% 
  group_by(ID_code, Plate_num) %>% 
  summarise(
    FOSA = sum(Conc_Corr_num, na.rm = TRUE),
    ID_code = first(ID_code),
    SciName = abbr_binom(first(SciName)),
    Coast = first(Coast), 
    Plate_num = first(Plate_num),
    .groups = 'drop'  # Suppress the grouping message
  ) %>% 
  mutate(
    Mn_or_not = ifelse(SciName == "M. novaeangliae", T, F)
  )

FOSA_baleen_idv_summ <- FOSA_baleen_idv %>% 
  group_by(Mn_or_not) %>% 
  summarise(
    Med_tot_PFAS = median(FOSA),
    SD_tot_PFAS = sd(FOSA)
  )

# Perform Kruskal-Wallis test
FOSA_by_ind <- kruskal.test(FOSA ~ Mn_or_not, data = FOSA_baleen_idv)

# Print the result
print(FOSA_by_ind)


