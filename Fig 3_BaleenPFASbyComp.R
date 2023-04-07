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
         Sample_seq ==1) %>% 
  group_by(ID_code, Compound) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup()



# Define the order of the compounds
PFASofInterest <- c("PFHxA", "PFHxS", "PFOA", "PFOS", "FOSA", "PFNA", 
                    "PFDA", "PFUdA", "PFDoA", "PFtrDA", "PFTeDA")

# Reorder the Compound variable according to the specified order
PFASBaleenCompSum$Compound <- factor(PFASBaleenCompSum$Compound, levels = PFASofInterest)

# Create the plot
Baleen_stacked_BP <- ggplot(PFASBaleenCompSum, 
                            aes(x = reorder(ID_code, -Total_PFAS), 
                                y = Total_PFAS, 
                                fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g in newest baleen)", fill = "Compound") +
  scale_fill_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #guides(fill = guide_legend(reverse = TRUE))
Baleen_stacked_BP

ggsave("Baleen_stacked_BP.pdf")



# Create the proportional plot
Baleen_stacked_BP_prop <- ggplot(PFASBaleenCompSum, 
                            aes(x = reorder(ID_code, -Total_PFAS), 
                                y = Total_PFAS/sum(Total_PFAS), 
                                fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual", y = "Proportion PFAS in newest baleen", fill = "Compound") +
  scale_fill_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill = guide_legend(reverse = TRUE))
Baleen_stacked_BP_prop

ggsave("Baleen_stacked_BP_prop.pdf")





