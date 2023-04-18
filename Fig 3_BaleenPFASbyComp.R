--------------------------------------------------------------------------------
# Baleen-PFAS by PFAS Compounds figure
# Matt Savoca
# Started on: 4/5/23
--------------------------------------------------------------------------------
  
  
source("Util.R")
library(pals)
library(ggpubr)


# 10 compounds that show up in at least 40% of all samples, plus PFOA
PFASofInterest <- Congener_screen$Compound[1:11]
PFASofInterest <- c("PFHxA", "PFHxS", "PFOA", "PFOS", "FOSA", "PFNA", 
                    "PFDA", "PFUdA", "PFDoA", "PFtrDA", "PFTeDA")



PFASBaleenCompSum <- Baleen_PFAS_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type == "baleen",
         if_else(ID_code == "IFAW13-158Mn", Sample_seq == 2, Sample_seq == 1)) %>% 
  group_by(ID_code, Compound, Coast, SciName) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup()





# Reorder the Compound variable according to the specified order
PFASBaleenCompSum$Compound <- factor(PFASBaleenCompSum$Compound, levels = PFASofInterest)


SppOrder = c("Eubalaena glacialis", "Balaenoptera musculus", "Balaenoptera borealis", 
             "Balaenoptera physalus", "Balaenoptera acutorostrata", "Megaptera novaeangliae")

PFASBaleenCompSum$SciName <- factor(PFASBaleenCompSum$SciName, levels = SppOrder)

PFASBaleenCompSum$ID_code <- factor(PFASBaleenCompSum$ID_code, levels = unique(PFASBaleenCompSum$ID_code[order(PFASBaleenCompSum$SciName)]))

ID_code_levels <- c("COA03-Mn", "IFAW17-274Mn", "COA15-0611Mn", "IFAW16-227Mn", "IFAW19-297Mn",
                    "COA1415-Mn", "IFAW17-317Mn", "IFAW19-287Mn", "IFAW13-158Mn", "IFAW20-009Mn",
                    "COA20-0804Ba", "COA14-0717Ba", "COA20-0808Ba", "COA16-06098Bb", "IFAW17-182Eg",
                    "GLMN1-3", "GLBP1-3", "GLBW1-61-9")


Baleen_stacked_BP <- ggplot(PFASBaleenCompSum, 
                            aes(x = fct_relevel(ID_code, ID_code_levels), 
                                y = Total_PFAS, 
                                fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g in newest baleen)", fill = "Compound") +
  scale_fill_manual(values = Compounds_pal) +
  facet_grid(.~Coast, scales = "free", space ="free") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

Baleen_stacked_BP

ggsave("Baleen_stacked_BP.pdf")



# Create the proportional plot
Baleen_stacked_BP_prop <- ggplot(PFASBaleenCompSum, 
                                aes(x = fct_relevel(ID_code, ID_code_levels), 
                                y = Total_PFAS/sum(Total_PFAS), 
                                fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual", y = "Proportion PFAS in newest baleen", fill = "Compound") +
  facet_grid(.~Coast, scales = "free", space ="free") +
  scale_fill_manual(values = Compounds_pal) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
Baleen_stacked_BP_prop

ggsave("Baleen_stacked_BP_prop.pdf")





