# Prepping PFAS data for analysis ----

source("Util.R")

View(Baleen_PFAS_data_comb)

# Test which compounds are present in >40% of all samples
Congener_screen <- Baleen_PFAS_data_comb %>% 
  group_by(Compound) %>% 
  summarise(Incl_40 = sum(Detect_binary)) %>% 
  arrange(-Incl_40)


Congener_conc <- Baleen_PFAS_data_comb %>% 
  group_by(Compound, Sample_type, ID_code) %>% 
  filter(Compound %in% c(Congener_screen$Compound[1:11])) %>% 
  summarise(Med_conc = median(Conc_Corr_num, na.rm = T)) %>% 
  arrange(-Med_conc)



Long_test_baleen <- Baleen_PFAS_data_comb %>% 
  filter(
    ##### Compound %in% Congener_screen$Compound %>% slice(1:10),
    Compound %in% c(Congener_screen$Compound[1:11]),
    #ID_code == "COA16-06098Bb",
    Sample_type == "baleen"
         )
  


# Split data by factor variable
my_data_split <- split(Long_test_baleen, Long_test_baleen$ID_code)

# Loop through the list of split data and create separate plots
plot_list <- purrr::map(my_data_split, function(df) {
  Long_test_plot <- ggplot(df, aes(Sample_seq, Conc_Corr_num)) +
    #geom_smooth() +
    geom_point() + 
    facet_wrap(Compound~., scales = "free_y") +
    labs(title = df$ID_code,
         x = "Baleen sample (gumline to tip)",
         y = "Concentration in baleen (ng/g)"
    )
  ggsave(paste0(unique(df$ID_code), ".pdf"), Long_test_plot)
})



Gum_test_plot <- ggplot(
  filter(Baleen_PFAS_data_comb, Sample_type == "gum",
         Compound %in% c(Congener_screen$Compound[1:11])), 
  aes(ID_code, Conc_Corr_num)) +
  #geom_smooth() +
  geom_point() + 
  facet_wrap(Compound~., scales = "free_y") +
  labs(title = ,
       x = "Individual ID",
       y = "Concentration in gum (ng/g)"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
Gum_test_plot


ggsave("PFAS gum tissue conc.pdf", Gum_test_plot)



Gum_stacked_BP <- 
  Baleen_PFAS_data_comb %>% 
  filter(Sample_type == "gum",
         Compound %in% c(Congener_screen$Compound[1:11])) %>%
  group_by(ID_code, Compound) %>% 
  mutate(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Total_PFAS)) %>% 
  ggplot(aes(x = reorder(ID_code, -Total_PFAS),
             y = Total_PFAS, fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g in gum)", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Gum_stacked_BP

ggsave("PFAS gum tissue conc_stackedbar.pdf", Gum_stacked_BP)



Gum_stacked_BP_prop <- Baleen_PFAS_data_comb %>% 
  filter(Sample_type == "gum",
         Compound %in% c(Congener_screen$Compound[1:11])) %>%
  group_by(ID_code, Compound) %>% 
  mutate(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Total_PFAS)) %>% 
  ggplot(aes(x = reorder(ID_code, -Total_PFAS),
             y = Total_PFAS/sum(Total_PFAS), fill = Compound)) +
  geom_col(position = "fill") +
  labs(x = "Individual", y = "Total PFAS proportion in gum by compound", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Gum_stacked_BP_prop

ggsave("PFAS gum tissue conc_stackedbar_prop.pdf", Gum_stacked_BP_prop)



--------------------------------------------------------------------------------
  # Prelim data cleaning and plotting for the wet tissues ----
--------------------------------------------------------------------------------
  
  source("Util.R")


library(pals)

# Generate the cols25() palette
cols <- cols25()
# Extract the first 12 colors
first_12_cols <- head(cols, 12)
# View the first 12 colors
first_12_cols


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



PFASWTCompSum <- PFAS_full_data_comb %>% 
  filter(Compound %in% PFASWTofInterest,
         Sample_type %in% c("gum", "liver", "skin", "blubber")
  ) %>% 
  group_by(ID_code, Compound, Sample_type) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup()



# Define the order of the compounds
# PFASofInterest <- c("PFHxA", "PFHxS", "PFOA", "PFOS", "FOSA", "PFNA", 
#                     "PFDA", "PFUdA", "PFDoA", "PFtrDA", "PFTeDA")

# Reorder the Compound variable according to the specified order
PFASBaleenCompSum$Compound <- factor(PFASBaleenCompSum$Compound, levels = PFASofInterest)

# Create the plot
WT_stacked_BP <- ggplot(PFASWTCompSum, 
                        aes(x = reorder(ID_code, -Total_PFAS), 
                            y = Total_PFAS, 
                            fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g)", fill = "Compound") +
  facet_wrap(.~Sample_type, scales = "free") +
  scale_fill_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
WT_stacked_BP





# Create the proportional plot

WT_stacked_BP_prop <- ggplot(PFASWTCompSum, 
                             aes(x = reorder(ID_code, -Total_PFAS), 
                                 y = Total_PFAS/sum(Total_PFAS), 
                                 fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual", y = "Total PFAS (ng/g)", fill = "Compound") +
  facet_wrap(.~Sample_type, scales = "free") +
  scale_fill_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
WT_stacked_BP_prop








#Junk code below here

Long_test_plot <- ggplot(Long_test_baleen, aes(Sample_seq, Conc_Corr_num)) +
  #geom_smooth() +
  geom_point() + 
  facet_wrap(Compound~., scales = "free_y") +
  labs(title = Long_test_baleen$ID_code,
       x = "Baleen sample (gumline to tip)",
       y = "Concentration in baleen (ng/g)"
  )


Long_test_plot





