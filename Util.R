# Utilities file for data cleaning and analysis

install.packages(c("tidyverse", "readxl", "NADA", "NADA2",))

library(tidyverse)
library(readxl)
library(NADA)
library(NADA2)
library(lme4)


Baleen_PFAS_samples <- read_xlsx("Baleen-PFAS master sample list.xlsx")

PFAS_data_raw <- read_csv("20220223_AreaCorrMassCorr.csv") %>% 
  rename("Sample_num" = "Sample_ID")


Baleen_PFAS_data_comb <- left_join(PFAS_data_raw, Baleen_PFAS_samples,
                                   by = c("Sample_num", "Exact_weight_g")) %>%
  mutate(Conc_Corr_num = as.numeric(ifelse(Conc_Corr == "<MDL", NA, Conc_Corr)),
         Censor_TF = as.logical(ifelse(Censored_Status == "Censored", "TRUE", "FALSE")),
         Detect_binary = ifelse(Conc_Corr_num %in% c(NA,0), 0, 1),
         Coast = ifelse(str_detect(ID_code, "GL"), "Pacific", "Atlantic"),
         Prey = ifelse(str_detect(ID_code, "(BW|Bb|BP|Eg)"), "Zooplankton", "Mixed (ZP/fish)")
  ) %>% 
  rename("Sample_seq" = "Sample_seq (if necessary)")



         
  
