#--------------------------------------------------------------------------------
  # Baleen-PFAS Total PFAS in baleen by 
  # Matt Savoca
  # Started on: 4/5/23
#--------------------------------------------------------------------------------


# Assuming you have defined shapes_custom and pal beforehand
  
# Assuming you have defined PFASofInterest previously
  
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
TL_PFAS_idv
  

TL_PFAS_idv_summ <- TL_PFAS_idv %>% 
  group_by(Prey) %>% 
  summarise(
    Med_tot_PFAS = median(Total_PFAS),
    SD_tot_PFAS = sd(Total_PFAS)
  )
    

# Perform Kruskal-Wallis test
kw_result <- kruskal.test(Total_PFAS ~ Prey, data = TL_PFAS_idv)

# Print the result
print(kw_result)
    
    
PFASbaleen_by_prey <- ggplot(TL_PFAS_idv, aes(x = Prey, y = Total_PFAS)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = SciName, shape = Coast), 
              width = 0.2, size = 4) +
  # scale_shape_manual(values = shapes_custom) +
  scale_color_manual(values = pal, 
                     guide = guide_legend(label.theme = element_text(angle = 0, face = "italic", size = 11)
                     )) +
  labs(y = bquote(paste("\u03A3"[10], "PFAS (ng/g d.w. in newest baleen)")),
       color = "Scientific name"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.45, 0.85),  # Adjust the position of the legend within the plot area (range 0 to 1)
        legend.justification = c(1, 1),  # Adjust the justification of the legend (range 0 to 1)
        legend.box.background = element_rect(color = "black"),  # Add a border (line) around the legend box
        legend.box.margin = margin(1, 1, 1, 1),
        axis.text.x = element_text(size = 14),  # Change the size of the x-axis tick text to 14
        axis.text.y = element_text(size = 14))  # Change the size of the y-axis tick text to 14
PFASbaleen_by_prey


ggsave("PFAS_by_prey_v2.jpg", PFASbaleen_by_prey, width = 6, height = 8)
  