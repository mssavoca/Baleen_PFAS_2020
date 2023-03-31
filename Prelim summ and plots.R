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





