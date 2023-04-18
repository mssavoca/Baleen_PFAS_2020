--------------------------------------------------------------------------------
  # Baleen-PFAS across tissue types 
  # Matt Savoca
  # Started on: 4/13/23
--------------------------------------------------------------------------------

source("Util.R")


PFASWTCompSum <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASWTofInterest,
         Sample_type %in% c("gum", "liver", "skin", "blubber")) %>% 
  group_by(ID_code, Compound, Sample_type) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup()



# Reorder the Compound variable according to the specified order
#PFASBaleenCompSum$Compound <- factor(PFASBaleenCompSum$Compound, levels = PFASofInterest)

# Create the plot
WT_stacked_BP <- ggplot(PFASWTCompSum, 
                        aes(x = reorder(ID_code, -Total_PFAS), 
                            y = Total_PFAS, 
                            fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g)", fill = "Compound") +
  facet_wrap(.~Sample_type, scales = "free") +
  scale_fill_manual(values = Compounds_pal) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
WT_stacked_BP

ggsave("WT_stacked_BP.pdf")



# Create the proportional plot

WT_stacked_BP_prop <- ggplot(PFASWTCompSum, 
                             aes(x = reorder(ID_code, -Total_PFAS), 
                                 y = Total_PFAS/sum(Total_PFAS), 
                                 fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual", y = "Proportion PFAS compounds", fill = "Compound") +
  facet_wrap(.~Sample_type, scales = "free") +
  scale_fill_manual(values = Compounds_pal) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
WT_stacked_BP_prop

ggsave("WT_stacked_BP_prop.pdf")

  
  

  
PFASCompSumAllTiss <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type != "gum-baleen interface",
         ifelse(Sample_type == "baleen", if_else(ID_code == "IFAW13-158Mn", 
                                                 Sample_seq == 2, Sample_seq == 1), TRUE)) %>% 
  group_by(ID_code, Sample_type) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  # mutate(PFAS_ratio = ifelse(Sample_type == "liver" & Compound %in% PFASofInterest,
  #                            Total_PFAS / sum(PFASCompSumAllTiss$Total_PFAS[PFASCompSumAllTiss$Sample_type == "baleen"]), 
  #                            NA)) %>% 
  ungroup()

  
  
PFAS_by_Tiss <- ggplot(filter(PFASCompSumAllTiss, 
                              ID_code %in% c("IFAW13-158Mn", "IFAW16-227Mn",
                                             "IFAW19-287Mn", "IFAW20-009Mn",
                                             "IFAW17-182Eg")), 
                          aes(x = reorder(Sample_type, -Total_PFAS), 
                              y = Total_PFAS, 
                              fill = Compound)) +
  geom_col() +
  labs(x = "Tissue type", y = "Total PFAS (ng/g)", fill = "Compound") +
  facet_wrap(.~ID_code, scales = "free") +
  scale_fill_manual(values = Compounds_pal) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PFAS_by_Tiss 

ggsave("PFAS_by_Tiss.pdf")



PFAS_by_Tiss_prop <- ggplot(filter(PFASCompSumAllTiss, 
                              ID_code %in% c("IFAW13-158Mn", "IFAW16-227Mn",
                                             "IFAW19-287Mn", "IFAW20-009Mn",
                                             "IFAW17-182Eg")), 
                       aes(x = reorder(Sample_type, -Total_PFAS), 
                           y = Total_PFAS/sum(Total_PFAS), 
                           fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Tissue type", y = "Proportion PFAS", fill = "Compound") +
  facet_wrap(.~ID_code, scales = "free") +
  scale_fill_manual(values = Compounds_pal) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PFAS_by_Tiss_prop 

ggsave("PFAS_by_Tiss_prop.pdf")










# Junk code below here 





Wet_tissues_congener_screen <- read_csv("20220409_Wet_DetFreq.csv")


PFAS_wet_tissues <- read_csv("20220409_Wet_AreaCorrMassCorr.csv") %>% 
  rename("Sample_num" = "Sample_Num") %>% 
  mutate(Sample_type = case_when(Matrix == "Skin" ~ "skin",
                                 Matrix == "Blubber" ~ "blubber",
                                 Matrix == "Liver" ~ "liver")) %>% 
  select(-Matrix, -`Sample Text`)



PFAS_wet_tissues_comb <- left_join(Baleen_PFAS_samples, PFAS_wet_tissues,
                                   by = c("Sample_num", "Exact_weight_g", "Sample_type")) %>%
  #filter(Sample_type %in% c("liver", "skin", "blubber")) %>% 
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



Wet_tissues_congener_screen <- PFAS_wet_tissues_comb %>% 
  group_by(Compound) %>% 
  summarise(Incl_40 = sum(Detect_binary),
            Overall_prop = Incl_40/21
  ) %>% 
  arrange(-Incl_40)




# 10 compounds that show up in at least 40% of all samples, plus PFOA
PFASWTofInterest <- Wet_tissues_congener_screen$Compound[1:9]
# PFASofInterest <- c("PFHxA", "PFHxS", "PFOA", "PFOS", "FOSA", "PFNA", 
#                     "PFDA", "PFUdA", "PFDoA", "PFtrDA", "PFTeDA")



