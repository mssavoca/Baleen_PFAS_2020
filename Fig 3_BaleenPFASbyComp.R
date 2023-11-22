--------------------------------------------------------------------------------
# Baleen-PFAS by PFAS Compounds figure
# Matt Savoca
# Started on: 4/5/23
--------------------------------------------------------------------------------
  
  
source("Util.R")
library(pals)
library(ggpubr)




PFASBaleenCompSum <- Baleen_PFAS_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type == "baleen",
         if_else(ID_code == "IFAW13-158Mn", Sample_seq == 2, Sample_seq == 1)) %>% 
  group_by(ID_code, ID_code_corr,Plate_num, Compound, Coast, SciName) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup() %>% 
  mutate(
    Combined_Column = ifelse(is.na(Plate_num), ID_code, 
                             paste(ID_code, ifelse(Plate_num == 2, "p2", "p1"), sep = " ")),
    Combined_Column_corr = ifelse(is.na(Plate_num), ID_code_corr, 
                                  paste(ID_code_corr, ifelse(Plate_num == 2, "p2", "p1"), sep = " "))
  )






# Reorder the Compound variable according to the specified order
PFASBaleenCompSum$Compound <- factor(PFASBaleenCompSum$Compound, levels = PFASofInterest)


SppOrder = c("Eubalaena glacialis", "Balaenoptera musculus", "Balaenoptera borealis", 
             "Balaenoptera physalus", "Balaenoptera acutorostrata", "Megaptera novaeangliae")

PFASBaleenCompSum$SciName <- factor(PFASBaleenCompSum$SciName, levels = SppOrder)

PFASBaleenCompSum$ID_code <- factor(PFASBaleenCompSum$ID_code, levels = unique(PFASBaleenCompSum$ID_code[order(PFASBaleenCompSum$SciName)]))

PFASBaleenCompSum$ID_code_corr <- factor(PFASBaleenCompSum$ID_code_corr, levels = unique(PFASBaleenCompSum$ID_code_corr[order(PFASBaleenCompSum$SciName)]))


Combined_Column_levels <- c("COA03-Mn p1", "IFAW17-274Mn p1", "COA15-0611Mn p1", "IFAW16-227Mn p1",
                            "IFAW19-287Mn p1", "IFAW19-287Mn p2", "COA1415-Mn p1", "IFAW17-317Mn p1", 
                            "IFAW13-158Mn p1", "IFAW20-009Mn p1", "COA20-0804Ba p1", "COA14-0717Ba p1", 
                            "COA20-0808Ba p1", "COA16-06098Bb p1", "IFAW17-182Eg p1",
                            "GLMN1-3 p1", "GLBP1-3 p1", "GLBW1-61-9 p1")
PFASBaleenCompSum$Combined_Column <- factor(PFASBaleenCompSum$Combined_Column, levels = Combined_Column_levels)



Combined_Column_levels_corr <- c("MH03602Mn p1", "IFAW17-274Mn p1", "COA150611Mn p1", "IFAW16-227Mn p1",
                            "IFAW19-287Mn p1", "IFAW19-287Mn p2", "COA141225Mn p1", "IFAW17-317Mn p1", 
                            "IFAW13-158Mn p1", "IFAW20-009Mn p1", "COA200804Ba p1", "COA140717Ba p1", 
                            "COA200808Ba p1", "COA150609Bb p1", "IFAW17-182Eg p1",
                            "GLMN1-3 p1", "GLBP1-3 p1", "GLBW1-61-9 p1")
PFASBaleenCompSum$Combined_Column_corr <- factor(PFASBaleenCompSum$Combined_Column_corr, levels = Combined_Column_levels_corr)


# Filter out "7:3 FTCA" compound
filtered_data <- PFASBaleenCompSum %>%
  filter(Compound != "7:3 FTCA")

# Create the stacked bar plot
Baleen_stacked_BP <- ggplot(filtered_data, 
                            aes(x = fct_relevel(Combined_Column_corr, Combined_Column_levels_corr), 
                                y = Total_PFAS, 
                                fill = Compound)) +
  geom_col() +
  labs(x = "", 
       y = bquote(paste("\u03A3"[10], "PFAS (ng/g d.w. in newest baleen)")),
       fill = "Compound") +
  scale_fill_manual(values = Compounds_pal) +
  facet_grid(.~Coast, scales = "free", space = "free") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank())

Baleen_stacked_BP

ggsave("Baleen_stacked_BP_v2.1.png", Baleen_stacked_BP, width = 12, height = 7)



# Create the proportional plot
Baleen_stacked_BP_prop <- ggplot(filtered_data, 
                                 aes(x = fct_relevel(Combined_Column_corr, Combined_Column_levels_corr), 
                                y = Total_PFAS/sum(Total_PFAS), 
                                fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual and plate", y = "Proportion PFAS in newest baleen", fill = "Compound") +
  facet_grid(.~Coast, scales = "free", space ="free") +
  scale_fill_manual(values = Compounds_pal) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
Baleen_stacked_BP_prop

ggsave("Baleen_stacked_BP_prop_v2.1.png", Baleen_stacked_BP_prop, width = 12, height = 8)










