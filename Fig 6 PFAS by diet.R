--------------------------------------------------------------------------------
  # Baleen-PFAS Total PFAS in baleen by 
  # Matt Savoca
  # Started on: 4/5/23
--------------------------------------------------------------------------------
  



# Assuming you have defined shapes_custom and pal beforehand
  
  TL_PFAS_idv <- Baleen_PFAS_data_comb %>%
    filter(
      Compound %in% c(PFASofInterest),
      Sample_seq == 1
    ) %>% 
    group_by(ID_code) %>% 
    summarise(
      Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
      ID_code = first(ID_code),
      SciName = abbr_binom(first(SciName)),
      Troph_level = first(Troph_level),
      Prey = first(Prey),
      Coast = first(Coast)
    ) %>%
    mutate(Prey = ifelse(Prey == "Mixed (ZP/fish)", "Zooplankton and fish",
                         ifelse(Prey == "Zooplankton", "Zooplankton-specialist", Prey))) %>%
    mutate(Prey = factor(Prey, levels = c("Zooplankton-specialist", "Zooplankton and fish"))) %>%
    ggplot(aes(x = Prey, y = Total_PFAS)) + 
    geom_boxplot() +
    geom_jitter(aes(color = SciName, shape = Coast), 
                width = 0.2, size = 4) +
    # scale_shape_manual(values = shapes_custom) +
    scale_color_manual(values = pal, 
                       guide = guide_legend(label.theme = element_text(angle = 0, face = "italic", size = 11)
                       )) +
    labs(y = bquote(paste("\u03A3"[11], "PFAS (ng/g d.w. in newest baleen)")),
         color = "Scientific name"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = c(0.45, 0.85),  # Adjust the position of the legend within the plot area (range 0 to 1)
          legend.justification = c(1, 1),  # Adjust the justification of the legend (range 0 to 1)
          legend.box.background = element_rect(color = "black"),  # Add a border (line) around the legend box
          legend.box.margin = margin(1, 1, 1, 1),
          axis.text.x = element_text(size = 14),  # Change the size of the x-axis tick text to 14
          axis.text.y = element_text(size = 14))  # Change the size of the y-axis tick text to 14
TL_PFAS_idv


ggsave("PFAS_by_prey.png", width = 6, height = 8)
  