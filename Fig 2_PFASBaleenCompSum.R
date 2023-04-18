--------------------------------------------------------------------------------
  # Baleen-PFAS full plate by summed PFAS figure
  # Matt Savoca
  # Started on: 4/5/23
--------------------------------------------------------------------------------
  
  
source("Util.R")
library(pals)
library(ggpubr)





# 10 comopounds that show up in at least 40% of all samples, plus PFOA
PFASofInterest <- Congener_screen$Compound[1:11]
PFASofInterest <- c("PFHxA", "PFHxS", "PFOA", "PFOS", "FOSA", "PFNA", "PFDA", "PFUdA", "PFDoA", "PFtrDA", "PFTeDA")



PFASBaleenCompSum <- Baleen_PFAS_data_comb %>% 
  filter(Compound %in% PFASofInterest,
         Sample_type == "baleen",
         !is.na(Sample_seq),
         is.finite(Sample_seq)) %>% 
  group_by(ID_code, Sample_seq) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup()

#View(PFASBaleenCompSum)

filter () %>% 
  select () %>% 
  summarise() %>% 
  mutate() %>% 
  arrang
# Plot without blue whale because only have one observation
# Plot without two Minke whales that have only one data point
PFAS_BCS_woBwMn <-  ggplot(filter(PFASBaleenCompSum, 
                                  !SciName %in% c("Balaenoptera musculus", "Megaptera novaeangliae"),
                                  !ID_code %in% c("COA20-0804Ba", "COA20-0808Ba")),
                           aes(x = -Sample_seq, y = Total_PFAS, 
                               color = abbr_binom(SciName),
                               shape = ID_code)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3.5) + 
  labs(x = "", 
       y = "Total PFAS (ng/g baleen)", 
       #color = "Species",
       shape = "Individual") +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  scale_shape_manual(values = shapes_custom) +
    scale_x_continuous(expand = c(0,1), 
                       labels = function(x) round(abs(x))
                       ) +
  facet_wrap(.~abbr_binom(SciName), scales = "free") +
  #scale_color_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 20) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "none")
PFAS_BCS_woBwMn

ggsave("PFAS_BCS_woBwMn.pdf")


# Plot without "IFAW13-158Mn" because has only have one observation
PFAS_BCS_MnOnly <- ggplot(filter(PFASBaleenCompSum, 
                                 SciName =="Megaptera novaeangliae",
                                 ID_code != "IFAW13-158Mn"),
                          aes(x = -Sample_seq, y = Total_PFAS, color = abbr_binom(SciName))) + 
  geom_line(size = 1.25) + 
  geom_point(size = 3) + 
  labs(x = "Baleen sample sequence", 
       y = "Total PFAS (ng/g baleen)") + 
  facet_wrap(.~ID_code, scales = "free", ncol = 5, nrow = 2) +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  scale_x_continuous(labels = function(x) round(abs(x)), limits = c(NA, -1), 
                     breaks = seq(-1, -100, -4)
                     ) + # Set x-axis limits and ticks
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

PFAS_BCS_MnOnly

ggsave("PFAS_BCS_MnOnly.pdf")


PFAS_BCS_comb <- ggarrange(PFAS_BCS_woBwMn, PFAS_BCS_MnOnly,
                           labels = c("A", "B"),
                           ncol = 1, nrow = 2,
                           padding = 0.1)
PFAS_BCS_comb
  