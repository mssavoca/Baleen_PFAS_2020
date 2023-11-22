# Utilities file for data cleaning and analysis

#install.packages(c("tidyverse", "readxl", "NADA", "NADA2",))

library(tidyverse)
library(readxl)
library(NADA)
library(NADA2)
library(lme4)
library(pals)
library(ggbreak)
library(ggforce)


# Generate the cols25() palette
cols <- cols25()
# Extract the first 12 colors
first_12_cols <- head(cols, 12)
# View the first 12 colors
first_12_cols


# color palette for figures
pal <- c("E. glacialis" = "#ff9f04", "B. acutorostrata" = "firebrick3", 
         "B. borealis" = "goldenrod2", "M. novaeangliae" = "gray40", 
         "B. physalus" = "chocolate3", "B. musculus" = "dodgerblue2")


shapes_custom <- c("IFAW13-158Mn" = 19, "IFAW20-009Mn" = 15,  "GLBW1-61-9" = 19,
                  "GLBP1-3" = 19, "IFAW19-287Mn" = 17, "IFAW17-182Eg" = 19,  
                  "IFAW16-227Mn" = 5, "IFAW17-274Mn" = 6, "IFAW17-317Mn" = 7,  
                  "GLMN1-3" = 8, "COA16-06098Bb" = 19, "COA15-0611Mn" = 9,  
                  "COA20-0808Ba" = 19, "COA20-0804Ba" = 15, "IFAW19-297Mn" = 10,
                  "COA14-0717Ba" = 17,  "COA1415-Mn" = 12,  "COA03-Mn" = 13)

Compounds_pal <- c( "PFHxA" = "#1F78C8", "PFHxS" = "#ff0000", "PFOA" = "#33a02c",
                  "PFOS" = "#6A33C2", "FOSA" = "#ff7f00", "PFNA" = "#565656",
                  "PFDA" = "#FFD700", "PFUdA" = "#a6cee3", "PFDoA" = "#FB6496", 
                  "PFtrDA" = "#b2df8a", "PFTeDA" = "#CAB2D6", "7:3 FTCA" = "#C8308C")

# Define the order of the compounds
PFASofInterest <- c("PFHxA", "PFOA", "PFOS", "FOSA", "PFNA", 
                    "PFDA", "PFUdA", "PFDoA", "PFtrDA", "PFTeDA", 
                    "7:3 FTCA")



# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}


# Import, clean and combine all data ----

# Sample metadata----
Baleen_PFAS_samples <- read_xlsx("Baleen-PFAS master sample list.xlsx")


# Wet tissue data----
Wet_tissues_congener_screen <- read_xlsx("20231010_Wet_DetFreq.xlsx")




PFAS_wet_tissues <- read_csv("20231010_Wet_AreaCorrMassCorr_ProbCmpdsRemoved.csv") %>% 
  rename("Sample_num" = "Sample_Num") %>% 
  # mutate(Sample_type = case_when(Matrix == "Skin" ~ "skin",
  #                                Matrix == "Blubber" ~ "blubber",
  #                                Matrix == "Liver" ~ "liver")) %>% 
  select(-Matrix)



PFAS_wet_tissues_comb <- left_join(Baleen_PFAS_samples, PFAS_wet_tissues,
                                   by = c("Sample_num", "Exact_weight_g")) %>%
  filter(Sample_type %in% c("liver", "skin", "blubber")) %>% 
  mutate(Conc_Corr_num = as.numeric(ifelse(Conc_Corr == "<MDL", NA, Conc_Corr)),
         Censor_TF = as.logical(ifelse(Censored_Status == "Censored", "TRUE", "FALSE")),
         Detect_binary = ifelse(Conc_Corr_num %in% c(NA,0), 0, 1),
         Coast = ifelse(str_detect(ID_code, "GL"), "Pacific", "Atlantic"),
         Prey = ifelse(str_detect(ID_code, "(BW|Bb|BP|Eg)"), "Zooplankton", "Mixed (ZP/fish)"), 
         SciName = case_when(Species == "Mn" ~ "Megaptera novaeangliae",
                             Species == "Ba" ~ "Balaenoptera acutorostrata",
                             Species == "Bb" ~ "Balaenoptera borealis",
                             Species == "Bw" ~ "Balaenoptera musculus",
                             Species == "Bp" ~ "Balaenoptera physalus",
                             Species == "Eg" ~ "Eubalaena glacialis")
  )






# Dry tissue (baleen and gum) data ---
PFAS_data_raw <- read_csv("20231010_AreaCorrMassCorr_ProbCmpdsRemoved.csv") %>% 
  rename("Sample_num" = "Sample_ID")


PFAS_data_raw <- read_csv("20231010_AreaCorrMassCorr_ProbCmpdsRemoved.csv") %>% 
  mutate(Sample_num = as.numeric(Sample_ID))


Baleen_PFAS_data_comb <- left_join(PFAS_data_raw, Baleen_PFAS_samples,
                                   by = c("Sample_num", "Exact_weight_g")) %>%
  mutate(Conc_Corr_num = as.numeric(ifelse(Conc_Corr == "<MDL", NA, Conc_Corr)),
         Censor_TF = as.logical(ifelse(Censored_Status == "Censored", "TRUE", "FALSE")),
         Detect_binary = ifelse(Conc_Corr_num %in% c(NA,0), 0, 1),
         Coast = ifelse(str_detect(ID_code, "GL"), "Pacific", "Atlantic"),
         Prey = ifelse(str_detect(ID_code, "(BW|Bb|BP|Eg)"), "Zooplankton", "Mixed (ZP/fish)"), 
         SciName = case_when(Species == "Mn" ~ "Megaptera novaeangliae",
                             Species == "Ba" ~ "Balaenoptera acutorostrata",
                             Species == "Bb" ~ "Balaenoptera borealis",
                             Species == "Bw" ~ "Balaenoptera musculus",
                             Species == "Bp" ~ "Balaenoptera physalus",
                             Species == "Eg" ~ "Eubalaena glacialis"),
         Troph_level = case_when(Species == "Mn" ~ 3.6,
                                Species == "Ba" ~ 3.4,
                                Species == "Bb" ~ 3.4,
                                Species == "Bw" ~ 3.2,
                                Species == "Bp" ~ 3.4,
                                Species == "Eg" ~ 3.2),
         ID_code_corr = ifelse(ID_code == "COA03-Mn", "MH03602Mn",
                               if_else(ID_code == "COA16-06098Bb", "COA150609Bb",
                                       ifelse(ID_code == "COA14-0717Ba", "COA140717Ba",
                                              ifelse(ID_code == "COA20-0804Ba", "COA200804Ba",
                                                     ifelse(ID_code == "COA20-0808Ba", "COA200808Ba",
                                                            ifelse(ID_code == "COA20-0804Ba", "COA200804Ba",
                                                                   ifelse(ID_code == "COA1415-Mn", "COA141225Mn",
                                                                          ifelse(ID_code == "COA15-0611Mn", "COA150611Mn", ID_code))))))))
                               ) %>% 
  rename("Sample_seq" = "Sample_seq (if necessary)") %>% 
  group_by(ID_code, Sample_seq) %>%
  mutate(Sample_seq_reversed = rev(seq_along(Sample_seq))) %>%  
  ungroup()

saveRDS(Baleen_PFAS_data_comb, file = "Baleen_PFAS_data_comb.RDS")



# Combine all cleaned datasets----
PFAS_full_data_comb <- merge(Baleen_PFAS_data_comb, PFAS_wet_tissues_comb, all = TRUE)





# Test which compounds are present in >40% of all samples
Dry_tissues_congener_screen <- Baleen_PFAS_data_comb %>% 
  filter(Sample_type != "gum") %>% 
  group_by(Compound) %>% 
  summarise(Incl_40 = sum(Detect_binary),
            Overall_prop = Incl_40/220) %>%  #just dry tissue samples here
  arrange(-Incl_40)
  

Wet_tissues_congener_screen <- PFAS_wet_tissues_comb %>% 
  group_by(Compound) %>% 
  summarise(Incl_40 = sum(Detect_binary),
            Overall_prop = Incl_40/21
  ) %>% 
  arrange(-Incl_40)



# 10 compounds that show up in at least 40% of all samples, plus PFOA
PFASWTofInterest <- Wet_tissues_congener_screen$Compound[1:8]
PFASWTofInterest <- c("PFNA", "PFOS", "PFUdA", "7:3 FTCA", "FOSA", "PFtrDA", 
                      "PFDA", "PFDoA", "PFHxA", "PFOA")


