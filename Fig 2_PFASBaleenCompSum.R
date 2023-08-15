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
  group_by(ID_code, Sample_seq, Plate_num) %>% 
  summarise(Total_PFAS = sum(Conc_Corr_num, na.rm = TRUE),
            SciName = first(SciName)) %>% 
  mutate(Plate_num = as.factor(Plate_num)) %>% 
  arrange(-Total_PFAS) %>% 
  ungroup()

#View(PFASBaleenCompSum)

filter () %>% 
  select () %>% 
  summarise() %>% 
  mutate() %>% 
  arrange
# Plot without blue whale because only have one observation
# Plot without two Minke whales that have only one data point
PFAS_BCS_woBwMn <-  ggplot(filter(PFASBaleenCompSum, 
                                  !SciName %in% c("Balaenoptera musculus", "Megaptera novaeangliae"),
                                  !ID_code %in% c("COA20-0804Ba", "COA20-0808Ba")),
                           aes(x = -Sample_seq, y = Total_PFAS, 
                               color = abbr_binom(SciName),
                               shape = ID_code)) + 
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3.5) + 
  labs(x = "", 
       y = bquote(paste("\u03A3"[11], "PFAS (ng/g d.w.)")),
       #color = "Species",
       shape = "Individual") +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  scale_shape_manual(values = shapes_custom) +
    scale_x_continuous(expand = c(0,1), 
                       labels = function(x) round(abs(x))
                       ) +
  facet_wrap(.~ID_code, scales = "free") +
  #scale_color_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 20) +
  theme(
    #strip.text = element_text(face = "italic"),
        legend.position = "none")
PFAS_BCS_woBwMn

ggsave("PFAS_BCS_woBwMn_v2.png", PFAS_BCS_woBwMn, width = 10, height = 7)


# Plot without "IFAW13-158Mn" because has only have one observation
PFAS_BCS_MnOnly <- ggplot(filter(PFASBaleenCompSum, 
                                 SciName =="Megaptera novaeangliae",
                                 !ID_code %in% c("IFAW13-158Mn", "IFAW16-227Mn",
                                                 "IFAW17-274Mn", "IFAW17-317Mn")),
                          aes(x = -Sample_seq, y = Total_PFAS, color = abbr_binom(SciName),
                              linetype = Plate_num, shape = Plate_num)) + 
  geom_line(size = 1) + 
  geom_point(size = 2.5) + 
  labs(x = "Baleen sample sequence", 
       y = bquote(paste("\u03A3"[11], "PFAS (ng/g d.w.)"))) +
  facet_wrap(.~ID_code, scales = "free", ncol = 3, nrow = 2) +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  scale_x_continuous(labels = function(x) round(abs(x)), limits = c(NA, -1), 
                     breaks = seq(-1, -100, -4)
                     ) + # Set x-axis limits and ticks
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

PFAS_BCS_MnOnly

ggsave("PFAS_BCS_MnOnly_v2.png", PFAS_BCS_MnOnly, width = 14, height = 6)


PFAS_BCS_comb <- ggarrange(PFAS_BCS_woBwMn, PFAS_BCS_MnOnly,
                           labels = c("A", "B"),
                           ncol = 1, nrow = 2,
                           padding = 0.1)
PFAS_BCS_comb






# PFAS:FOSA ratio plots along the plate ----


# Function to create ggplot for each ratio column and save the plot
create_ratio_plot <- function(data, ratio_column) {
  plot <- ggplot(data, aes(x = -Sample_seq, y = .data[[ratio_column]])) +
    geom_point() +
    geom_line() +
    labs(x = "Sample_seq", y = ratio_column) +
    facet_wrap(Plate_num ~ ID_code, scales = "free") +
    scale_x_continuous(labels = function(x) round(abs(x)), limits = c(NA, -1), 
                       breaks = seq(-1, -100, -4))
  
  # Generate a unique file name based on the ratio column
  file_name <- paste0(ratio_column, ".png")
  
  # Save the plot as a PNG file
  ggsave(plot, file = file_name, width = 10, height = 7.5, dpi = 300)
  
  # Return the plot
  return(plot)
}

# Define the ratio columns
ratio_columns <- c("Ratio_PFOS_FOSA", "Ratio_PFOS_PFUdA", "Ratio_FOSA_PFUdA")

saveRDS(Baleen_PFAS_ratios, "Baleen_PFAS_ratios.rds")


# Filter the data for ID_code with more than three Sample_seq values
filtered_data <- Baleen_PFAS_ratios %>%
  group_by(ID_code) %>%
  filter(
    !ID_code %in% c("COA16-06098Bb", "GLBP1-3"), # removing individuals with lots of zeros across plots
    n_distinct(Sample_seq) > 3) %>%
  ungroup()

# Apply the map function to create plots for each ratio column and save the plots
plots <- map(ratio_columns, ~create_ratio_plot(filtered_data, .x))

# Display the plots
plots





# Zoomed in plot on specific whales of interest----

#NARW calf baleen plot


PFAS_NARW_only <-  ggplot(filter(PFASBaleenCompSum, 
                                 ID_code == "IFAW17-182Eg"),
                           aes(x = -Sample_seq, y = Total_PFAS, 
                               color = abbr_binom(SciName),
                               shape = ID_code)) + 
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3.5) + 
  labs(x = "Baleen sample sequence", 
       y = "Total PFAS (ng/g baleen)", 
       #color = "Species",
       shape = "Individual") +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  scale_shape_manual(values = shapes_custom) +
  scale_x_continuous(expand = c(0,1), 
                     labels = function(x) round(abs(x))
  ) +
  labs(y = bquote(paste("\u03A3"[11], "PFAS (ng/g d.w.)"))) +
  facet_wrap(.~ID_code, scales = "free") +
  #scale_color_manual(values=as.vector(cols25(25))) +
  theme_bw(base_size = 20) +
  theme(
    #strip.text = element_text(face = "italic"),
        legend.position = "none")
PFAS_NARW_only

ggsave("PFAS_NARW_only.png", PFAS_NARW_only, width = 9, height = 5.5)



#Spinnaker's baleen plot

Spinnaker_df <- PFASBaleenCompSum %>% 
  filter(SciName == "Megaptera novaeangliae",
         ID_code == "COA15-0611Mn",
         Sample_seq < 11 
         ) %>% 
  select(-Plate_num) %>% 
  left_join(filtered_data, by = c("ID_code", "Sample_seq"))



PFAS_Spinnaker_baleen <- ggplot(Spinnaker_df,
                          aes()) + 
  geom_line(aes(x = -Sample_seq, y = Total_PFAS, color = abbr_binom(SciName)),
            size = 1) + 
  geom_point(aes(x = -Sample_seq, y = Total_PFAS, color = abbr_binom(SciName)),
             size = 3) +
  labs(x = "", 
       y = bquote(paste("\u03A3"[11], "PFAS (ng/g d.w.)"))
       ) + 
  facet_wrap(.~ID_code, scales = "free") +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  scale_x_continuous(labels = function(x) round(abs(x)), limits = c(-10.5, 0), 
                     breaks = seq(0, -100, -4)
  ) + # Set x-axis limits and ticks
  ylim(10,40) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

PFAS_Spinnaker_baleen

ggsave("PFAS_Spinnaker_baleen.png", PFAS_Spinnaker_baleen, width = 9, height = 5)



original_vector <- 0:20

# Multiply each value by 1.25
result_vector <- original_vector * 1.25

# Append the original vector to the result vector
# Combine vectors as separate columns
combined_matrix <- cbind(result_vector, original_vector)

# Convert the matrix to a data frame
combined_dataframe <- as.data.frame(combined_matrix)

# Rename columns
colnames(combined_dataframe) <- c("Result", "Original")

# Print the combined data frame
print(combined_dataframe)


PFAS_Spinnaker_FOSA_PFOS_ratio <- ggplot(Spinnaker_df,
                                aes()) + 
  geom_line(aes(x = -Sample_seq, y = 1/Ratio_PFOS_FOSA),
            linetype = "dashed", size = 1, color = "gray40") + 
  geom_point(aes(x = -Sample_seq, y = 1/Ratio_PFOS_FOSA),
             shape = 17, size = 4, color = "gray40") +
  labs(x = "Baleen sample sequence", 
       y = "FOSA:PFOS ratio") + 
  scale_x_continuous(labels = function(x) round(abs(x)), limits = c(-10.5, 0), 
                     breaks = seq(0, -100, -4)
  ) + # Set x-axis limits and ticks
  ylim(0,3) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")
PFAS_Spinnaker_FOSA_PFOS_ratio


ggsave("PFAS_Spinnaker_FOSA_PFOS_ratio.png", PFAS_Spinnaker_FOSA_PFOS_ratio, width = 8.4, height = 4.5)




PFAS_Spinnaker_comb <- ggarrange(PFAS_Spinnaker_baleen, PFAS_Spinnaker_FOSA_PFOS_ratio,
                           labels = c("A", "B"),
                           ncol = 1, nrow = 2,
                           padding = 0.1)
PFAS_Spinnaker_comb






  