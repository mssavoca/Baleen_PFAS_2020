--------------------------------------------------------------------------------
# PFAS prelim data analysis
# Matt Savoca
# Started on: 3/15/23
--------------------------------------------------------------------------------

source("Util.R")

## PFAS congeners to include in analysis
Congener_screen <- Baleen_PFAS_data_comb %>% 
  group_by(Compound) %>% 
  summarise(Incl_40 = sum(Detect_binary)) %>% 
  arrange(-Incl_40)




View(Baleen_PFAS_data_comb)




# anova for PFAS gums different krill feeding whales vs fish-feeding whales

# linear regressions within individual plate use RAW data (not NADA-corrected)
# look into log-transformations bc data is likely not normally distributed (ie use non-parametric test)

# gum PFAS burden vs baleen PFAS burden (NADA-corrected);

# for single point comparisons use ***CORRECTED DATA*** no <MDL values can be used 
#1) compare ratio most recent baleen vs. gum
# check ratio of baleen vs. gum over time on the plate


#for me: correct wet tissue data; analyze dry tissues


# check ratio of baleen vs. gum over time on the plate
  
  
Baleen_vs_gum <-  Baleen_PFAS_data_comb %>%
    filter(Sample_type == "gum",
           Compound == "PFOS") %>% 
    drop_na(Conc_Corr_num) %>% 
    group_by(ID_code, Compound) %>%
  summarise(G_B_Ratio = Conc_Corr_num[Sample_type == "gum"] / Conc_Corr_num[Sample_type == "baleen"],
            Ratio_text = Description[Sample_type == "baleen"])


# Code from MFC on 3/22/23

#gum to baleen ratio plot code----

get_gb_ratio <- function(sample, conc) {
  conc[sample == "gum"] / conc[sample == "baleen"]
}

valid_individuals <- Baleen_PFAS_data_comb %>% 
  drop_na(Sample_type, Conc_Corr_num) %>% 
  group_by(ID_code, Compound) %>% 
  summarize(is_valid = all(c("gum", "baleen") %in% Sample_type),
            .groups = "drop") %>% 
  filter(is_valid)

Gum_baleen_ratio_plot <- Baleen_PFAS_data_comb %>% 
  drop_na(Sample_type, Conc_Corr_num) %>% 
  semi_join(valid_individuals, by = c("ID_code", "Compound")) %>% 
  filter(Sample_type %in% c("gum", "baleen"),
         Compound %in% c(Congener_screen$Compound[1:11])) %>% 
  group_by(ID_code, Compound) %>% 
  summarize(bg_ratio = get_bg_ratio(Sample_type, Conc_Corr_num),
            baleen_seq = Sample_seq[Sample_type == "baleen"],
            .groups = "drop") %>% 
  ggplot(aes(baleen_seq, log10(bg_ratio), color = Compound)) +
  geom_line() +
  #geom_smooth(method = "lm", formula = y ~ x, se = FALSE, alpha = 0.5) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ ID_code, scales = "free_x") +
  scale_x_continuous(breaks = seq(1, 20, 2)) +
  labs(x = "Baleen sample (gumline to tip)",
       y = "Log10(gum : baleen ratio)"
  ) +
  theme_bw(base_size = 14)
Gum_baleen_ratio_plot
             
ggsave("Gum_baleen_ratio_plot.pdf")    




Baleen_long_plot <- Baleen_PFAS_data_comb %>% 
  drop_na(Sample_type, Conc_Corr_num, ID_code) %>% 
  filter(Sample_type == "baleen",
         Compound %in% c(Congener_screen$Compound[1:11])) %>% 
  group_by(ID_code, Compound) %>% 
  ggplot(aes(Sample_seq, log10(Conc_Corr_num), color = Compound)) +
  geom_line() +
  #geom_smooth(method = "lm", formula = y ~ x, se = FALSE, alpha = 0.5) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ ID_code, scales = "free_x") +
  scale_x_continuous(breaks = seq(1, 30, 2)) +
  labs(x = "Baleen sample (gumline to tip)",
       y = "Log10(concentration in baleen (ng/g))"
  ) +
  theme_bw(base_size = 16)
Baleen_long_plot

ggsave("Gum_baleen_ratio_plot.pdf")  



# NADA-corrected baleen analyses ----


# Calculate NADA-corrected data summaries

# vector with compounds of interest
PFASofinterest = Congener_screen$Compound[1:11]

# create a function to compute the cfit and return a data frame with the results
compute_cfit_df <- function(compound_name) {
  # subset the data frame for the current compound
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name) %>% 
    group_by(Compound)
  
  # compute the cfit for the current compound
  mycenfit <- cfit(compound_data$Conc_NoCorr + 0.0000001, compound_data$Censor_TF)
  
  # create a data frame with the compound name and the cfit coefficients
  result_df <- data.frame(Compound = compound_name,
                          Q10 = mycenfit$Q10,
                          Q25 = mycenfit$Q25,
                          Q50 = mycenfit$Q50,
                          Q75 = mycenfit$Q75,
                          Q90 = mycenfit$Q90,
                          stringsAsFactors = FALSE)
  
  # return the data frame
  return(result_df)
}

# use the map_df function to apply the compute_cfit_df function to each element of the PFASofinterest list and combine the results into a single data frame
mycenfit_df <- map_df(PFASofinterest, compute_cfit_df)

# convert character columns to numeric, replacing non-numeric values with NA
mycenfit_df2 <- mycenfit_df %>% 
  mutate_all(type.convert)

# convert Q10, Q25, and Q50 columns to numeric type and replace non-numeric values with NA
mycenfit_df2$Q10 <- ifelse(is.na(as.numeric(mycenfit_df2$Q10)), NA, as.numeric(mycenfit_df2$Q10))
mycenfit_df2$Q25 <- ifelse(is.na(as.numeric(mycenfit_df2$Q25)), NA, as.numeric(mycenfit_df2$Q25))
mycenfit_df2$Q50 <- ifelse(is.na(as.numeric(mycenfit_df2$Q50)), NA, as.numeric(mycenfit_df2$Q50))

# check the data types of each column again
str(mycenfit_df2)

# replace all NAs with 0
mycenfit_df2[is.na(mycenfit_df2)] <- 0


# Convert data frame to tidy long format
myquantiles_long <- gather(mycenfit_df2, key = "Quantile", value = "Value", -Compound, convert = TRUE)
myquantiles_long$Value = as.numeric(myquantiles_long$Value)

# Create box plot

Baleen_PFAS_grouped <-  ggplot(myquantiles_long, aes(x = fct_reorder(Compound, -Value, median), y = Value)) +
  geom_boxplot(aes(color = Compound), outlier.shape = NA) +
  geom_linerange(aes(ymin = Value, ymax = Value), size = 1) +
  labs(x = "Compound", y = "Concentration in baleen (ng/g)") +
  #scale_color_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
  guides(color = FALSE) +
  theme_classic(base_size = 14)
Baleen_PFAS_grouped   


ggsave("Baleen_PFAS_grouped.pdf")  
  



# Baleen-PFAS data by Compound and Individual----

# vector with compounds of interest
PFASofinterest = Congener_screen$Compound[1:11]


# Baleen_PFAS_data_comb_clip <- Baleen_PFAS_data_comb %>% 
#   filter(ID_code %in% c("IFAW20-009Mn", "COA03-Mn"))

# create a function to compute the cfit and return a data frame with the results
compute_cfit_df <- function(compound_name, id_code) {
  # subset the data frame for the current compound and ID_code
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name & ID_code == id_code)
  
  # check if there are at least two non-missing, uncensored, distinct values
  if (nrow(compound_data %>% filter(!is.na(Conc_NoCorr), Censor_TF == 0) %>% distinct()) < 2) {
    # return a data frame with NA values for the cfit coefficients
    result_df <- data.frame(ID_code = id_code,
                            Compound = compound_name,
                            Q10 = NA_real_,
                            Q25 = NA_real_,
                            Q50 = NA_real_,
                            Q75 = NA_real_,
                            Q90 = NA_real_,
                            stringsAsFactors = FALSE)
  } else {
    # compute the cfit for the current compound and ID_code
    mycenfit <- cfit(compound_data$Conc_NoCorr + 0.0000001, compound_data$Censor_TF)
    
    # create a data frame with the compound name, ID_code, and the cfit coefficients
    result_df <- data.frame(ID_code = id_code,
                            Compound = compound_name,
                            Q10 = as.numeric(mycenfit$Q10),
                            Q25 = as.numeric(mycenfit$Q25),
                            Q50 = as.numeric(mycenfit$Q50),
                            Q75 = as.numeric(mycenfit$Q75),
                            Q90 = as.numeric(mycenfit$Q90),
                            stringsAsFactors = FALSE)
  }
  
  # return the data frame
  return(result_df)
}



# create an empty data frame to store the results
mycenfit_df_ID <- data.frame(ID_code = character(), Compound = character(), Q10 = double(),
                          Q25 = double(), Q50 = double(), Q75 = double(), Q90 = double(),
                          stringsAsFactors = FALSE)

# loop over each ID_code and apply the compute_cfit_df function
for (id_code in unique(Baleen_PFAS_data_comb$ID_code)) {
  
  # filter data for the current ID_code and check if the Sample_seq is greater than 2
  if (Baleen_PFAS_data_comb %>%
      filter(ID_code == id_code, ) %>%
      distinct(Sample_seq) %>%
      nrow() > 2) {
    
    # use the map_df function to apply the compute_cfit_df function to each element of the PFASofinterest list for the current ID_code and combine the results into a single data frame
    mycenfit_df_id <- map_df(PFASofinterest, ~compute_cfit_df(.x, id_code))
    
    # append the results to the mycenfit_df data frame
    mycenfit_df_ID <- rbind(mycenfit_df_ID, mycenfit_df_id)
  } else {
    # print a message indicating that there is not enough data for the current ID_code
    message(paste0("Not enough data for ID_code ", id_code))
  }
}


# Save the mycenfit_df_ID object to a file
saveRDS(mycenfit_df_ID, "mycenfit_df_ID.rds")
mycenfit_df_ID <- readRDS("mycenfit_df_ID.rds")




# replace all NAs with 0
mycenfit_df_ID[is.na(mycenfit_df_ID)] <- 0


# Convert data frame to tidy long format
myquantiles_long_ID <- gather(mycenfit_df_ID, key = "Quantile", value = "Value", -c(Compound, ID_code), convert = TRUE)
myquantiles_long_ID$Value = as.numeric(myquantiles_long_ID$Value)

# Create box plot

Baleen_PFAS_grouped_ID <- ggplot(myquantiles_long_ID, aes(x = fct_reorder(Compound, -Value, median), y = Value)) +
  geom_boxplot(aes(color = Compound), outlier.shape = NA) +
  geom_linerange(aes(ymin = Value, ymax = Value), size = 1) +
  labs(x = "Compound", y = "Concentration in baleen (ng/g)") +
  facet_wrap(.~ID_code, scales = "free_y") +
  #scale_color_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
  guides(color = FALSE) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Baleen_PFAS_grouped_ID  

ggsave("Baleen_PFAS_grouped_ID.pdf")  

sum_baleen_vals <- mycenfit_df_ID %>% 
  group_by(ID_code) %>% 
  summarise(sum(Q50, na.rm = T))


Baleen_stacked_BP <- ggplot(filter(mycenfit_df_ID, ID_code != "IFAW17-274Mn"), 
                                   aes(x = reorder(ID_code, -Q50),
             y = Q50, fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g in gum)", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Baleen_stacked_BP

ggsave("Baleen_stacked_BP.pdf") 



Baleen_stacked_BP <- ggplot(filter(mycenfit_df_ID, 
                                   !ID_code %in% c("IFAW20-009Mn")), 
                            aes(x = reorder(ID_code, -Q50),
                                y = Q50, fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g in baleen)", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Baleen_stacked_BP

ggsave("Baleen_stacked_BP.pdf")  



Baleen_stacked_BP_prop <- ggplot(filter(mycenfit_df_ID, 
                                         !ID_code %in% c("IFAW20-009Mn")), 
                            aes(x = reorder(ID_code, -Q50),
                                y = Q50/sum(Q50), fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual", y = "Total PFAS (ng/g in baleen)", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Baleen_stacked_BP_prop

ggsave("Baleen_stacked_BP_prop.pdf") 





# Baleen-PFAS data by Compound and Coast----

# vector with compounds of interest
PFASofinterest = Congener_screen$Compound[1:11]


# Baleen_PFAS_data_comb_clip <- Baleen_PFAS_data_comb %>% 
#   filter(ID_code %in% c("IFAW20-009Mn", "COA03-Mn"))

# create a function to compute the cfit and return a data frame with the results
compute_cfit_df <- function(compound_name, coast) {
  # subset the data frame for the current compound and ID_code
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name & Coast == coast)
  
  # check if there are at least two non-missing, uncensored, distinct values
  if (nrow(compound_data %>% filter(!is.na(Conc_NoCorr), Censor_TF == 0) %>% distinct()) < 2) {
    # return a data frame with NA values for the cfit coefficients
    result_df <- data.frame(Coast == coast,
                            Compound = compound_name,
                            Q10 = NA_real_,
                            Q25 = NA_real_,
                            Q50 = NA_real_,
                            Q75 = NA_real_,
                            Q90 = NA_real_,
                            stringsAsFactors = FALSE)
  } else {
    # compute the cfit for the current compound and ID_code
    mycenfit <- cfit(compound_data$Conc_NoCorr + 0.0000001, compound_data$Censor_TF)
    
    # create a data frame with the compound name, ID_code, and the cfit coefficients
    result_df <- data.frame(Coast = coast,
                            Compound = compound_name,
                            Q10 = as.numeric(mycenfit$Q10),
                            Q25 = as.numeric(mycenfit$Q25),
                            Q50 = as.numeric(mycenfit$Q50),
                            Q75 = as.numeric(mycenfit$Q75),
                            Q90 = as.numeric(mycenfit$Q90),
                            stringsAsFactors = FALSE)
  }
  
  # return the data frame
  return(result_df)
}



# create an empty data frame to store the results
mycenfit_df_Coast <- data.frame(Coast = character(), Compound = character(), Q10 = double(),
                             Q25 = double(), Q50 = double(), Q75 = double(), Q90 = double(),
                             stringsAsFactors = FALSE)

# loop over each ID_code and apply the compute_cfit_df function
for (coast in unique(Baleen_PFAS_data_comb$Coast)) {
  
  # filter data for the current ID_code and check if the Sample_seq is greater than 2
  if (Baleen_PFAS_data_comb %>%
      filter(Coast == coast, ) %>%
      distinct(Sample_seq) %>%
      nrow() > 2) {
    
    # use the map_df function to apply the compute_cfit_df function to each element of the PFASofinterest list for the current ID_code and combine the results into a single data frame
    mycenfit_df_coast <- map_df(PFASofinterest, ~compute_cfit_df(.x, coast))
    
    # append the results to the mycenfit_df data frame
    mycenfit_df_Coast <- rbind(mycenfit_df_Coast, mycenfit_df_coast)
  } else {
    # print a message indicating that there is not enough data for the current ID_code
    message(paste0("Not enough data for Coast ", coast))
  }
}


# Save the mycenfit_df_ID object to a file
saveRDS(mycenfit_df_Coast, "mycenfit_df_Coast.rds")
mycenfit_df_Coast <- readRDS("mycenfit_df_Coast.rds")




# replace all NAs with 0
mycenfit_df_Coast[is.na(mycenfit_df_Coast)] <- 0


# Convert data frame to tidy long format
myquantiles_long_Coast <- gather(mycenfit_df_Coast, key = "Quantile", value = "Value", -c(Compound, Coast), convert = TRUE)
myquantiles_long_Coast$Value = as.numeric(myquantiles_long_Coast$Value)

# Create box plot

Baleen_PFAS_grouped_Coast <- ggplot(myquantiles_long_Coast, aes(x = fct_reorder(Compound, -Value, median), y = Value)) +
  geom_boxplot(aes(color = Compound), outlier.shape = NA) +
  geom_linerange(aes(ymin = Value, ymax = Value), size = 1) +
  labs(x = "Compound", y = "Concentration in baleen (ng/g)") +
  facet_wrap(.~Coast) +
  #scale_color_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
  guides(color = FALSE) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Baleen_PFAS_grouped_Coast  

ggsave("Baleen_PFAS_grouped_Coast.pdf") 




Baleen_stacked_BP_Coast <- ggplot(mycenfit_df_Coast, 
                            aes(x = reorder(Coast, -Q50),
                                y = Q50, fill = Compound)) +
  geom_col() +
  labs(x = "Individual", y = "Total PFAS (ng/g in gum)", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Baleen_stacked_BP_Coast

ggsave("Baleen_stacked_BP_Coast.pdf")  



Baleen_stacked_BP_prop <- ggplot(filter(mycenfit_df_ID, 
                                        !ID_code %in% c("IFAW17-274Mn", "IFAW20-009Mn")), 
                                 aes(x = reorder(ID_code, -Q50),
                                     y = Q50/sum(Q50), fill = Compound)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Individual", y = "Total PFAS (ng/g in gum)", fill = "Compound") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Baleen_stacked_BP_prop

ggsave("Baleen_stacked_BP_prop.pdf") 





 # Baleen-PFAS data by Compound and Prey----

# vector with compounds of interest
PFASofinterest = Congener_screen$Compound[1:11]


# Baleen_PFAS_data_comb_clip <- Baleen_PFAS_data_comb %>% 
#   filter(ID_code %in% c("IFAW20-009Mn", "COA03-Mn"))

# create a function to compute the cfit and return a data frame with the results
compute_cfit_df <- function(compound_name, prey) {
  # subset the data frame for the current compound and ID_code
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name & Prey == prey)
  
  # check if there are at least two non-missing, uncensored, distinct values
  if (nrow(compound_data %>% filter(!is.na(Conc_NoCorr), Censor_TF == 0) %>% distinct()) < 2) {
    # return a data frame with NA values for the cfit coefficients
    result_df <- data.frame(Prey == prey,
                            Compound = compound_name,
                            Q10 = NA_real_,
                            Q25 = NA_real_,
                            Q50 = NA_real_,
                            Q75 = NA_real_,
                            Q90 = NA_real_,
                            stringsAsFactors = FALSE)
  } else {
    # compute the cfit for the current compound and ID_code
    mycenfit <- cfit(compound_data$Conc_NoCorr + 0.0000001, compound_data$Censor_TF)
    
    # create a data frame with the compound name, ID_code, and the cfit coefficients
    result_df <- data.frame(Prey = prey,
                            Compound = compound_name,
                            Q10 = as.numeric(mycenfit$Q10),
                            Q25 = as.numeric(mycenfit$Q25),
                            Q50 = as.numeric(mycenfit$Q50),
                            Q75 = as.numeric(mycenfit$Q75),
                            Q90 = as.numeric(mycenfit$Q90),
                            stringsAsFactors = FALSE)
  }
  
  # return the data frame
  return(result_df)
}



# create an empty data frame to store the results
mycenfit_df_Prey <- data.frame(Prey = character(), Compound = character(), Q10 = double(),
                               Q25 = double(), Q50 = double(), Q75 = double(), Q90 = double(),
                               stringsAsFactors = FALSE)

# loop over each ID_code and apply the compute_cfit_df function
for (prey in unique(Baleen_PFAS_data_comb$Prey)) {
  
  # filter data for the current ID_code and check if the Sample_seq is greater than 2
  if (Baleen_PFAS_data_comb %>%
      filter(Prey == prey, ) %>%
      distinct(Sample_seq) %>%
      nrow() > 2) {
    
    # use the map_df function to apply the compute_cfit_df function to each element of the PFASofinterest list for the current ID_code and combine the results into a single data frame
    mycenfit_df_prey <- map_df(PFASofinterest, ~compute_cfit_df(.x, prey))
    
    # append the results to the mycenfit_df data frame
    mycenfit_df_Prey <- rbind(mycenfit_df_Prey, mycenfit_df_prey)
  } else {
    # print a message indicating that there is not enough data for the current ID_code
    message(paste0("Not enough data for Prey ", Prey))
  }
}


# Save the mycenfit_df_ID object to a file
saveRDS(mycenfit_df_Prey, "mycenfit_df_Prey.rds")
mycenfit_df_Prey <- readRDS("mycenfit_df_Prey.rds")




# replace all NAs with 0
mycenfit_df_Coast[is.na(mycenfit_df_Prey)] <- 0


# Convert data frame to tidy long format
myquantiles_long_Prey <- gather(mycenfit_df_Prey, key = "Quantile", value = "Value", -c(Compound, Prey), convert = TRUE)
myquantiles_long_Prey$Value = as.numeric(myquantiles_long_Prey$Value)

# Create box plot

Baleen_PFAS_grouped_Prey <- ggplot(myquantiles_long_Prey, aes(x = fct_reorder(Compound, -Value, median), y = Value)) +
  geom_boxplot(aes(color = Compound), outlier.shape = NA) +
  geom_linerange(aes(ymin = Value, ymax = Value), size = 1) +
  labs(x = "Compound", y = "Concentration in baleen (ng/g)") +
  facet_wrap(~fct_relevel(Prey, "Zooplankton", "Mixed (ZP/fish)")) +
  #scale_color_manual(values = c("#377eb8", "#4daf4a", "#984ea3", "#e41a1c", "#ff7f00")) +
  guides(color = FALSE) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Baleen_PFAS_grouped_Prey  

ggsave("Baleen_PFAS_grouped_Prey.pdf") 














#Junk code below here ----


RecentBaleen_vs_gum <- Baleen_PFAS_data_comb %>% 
  filter(ID_code %in% c("IFAW20-009Mn", "IFAW19-287Mn",
                        "IFAW17-182Eg", "IFAW16-227Mn", "IFAW17-274Mn",
                        "IFAW17-317Mn", "COA20-0804Ba"))
group_by(ID_code)









# Baleen-PFAS data by compound----
# create an empty list to store the results
mycenfit_list <- list()

# loop over the compound names and compute the cfit for each one
for (compound_name in PFASofinterest) {
  # subset the data frame for the current compound
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name)
  
  # compute the cfit for the current compound
  mycenfit <- cfit(compound_data$Conc_NoCorr + 0.0000001, compound_data$Censor_TF)
  
  # add the cfit result to the list
  mycenfit_list[[compound_name]] <- mycenfit
}


# create an empty list to store the results
myquantiles_list <- list()

# loop over the compound names and compute the quantiles for each one
for (compound_name in PFASofinterest) {
  # subset the data frame for the current compound
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name)
  
  # compute the cfit for the current compound
  mycenfit <- cfit(compound_data$Conc_NoCorr + 0.0000001, compound_data$Censor_TF)
  
  # extract the desired quantile values
  myquantiles <- c(Q10 = as.numeric(mycenfit$Q10), 
                   Q25 = as.numeric(mycenfit$Q25), 
                   Q50 = as.numeric(mycenfit$Q50), 
                   Q75 = as.numeric(mycenfit$Q75), 
                   Q90 = as.numeric(mycenfit$Q90))
  
  # add the quantile values to the list
  myquantiles_list[[compound_name]] <- myquantiles
}




# Define a new function to compute the cfit with missing and censored values excluded
cfit_mod <- function(y, censor, conf=0.95) {
  # subset the data to exclude missing and censored values
  ysub <- subset(y, complete.cases(y, censor))
  censor_sub <- subset(censor, complete.cases(y, censor))
  
  # compute the cfit for the subsetted data
  if (length(unique(ysub)) < 2) {
    return(list(point.est=NA, ci.lwr=NA, ci.upr=NA))
  } else {
    fit <- enparCensored(ysub, censor_sub, ci=TRUE, pivot.statistic="t", ci.sample.size=length(ysub))
    return(fit)
  }
}



# create an empty list to store the results
mycenfit_list_ID <- list()

# loop over the compound names and compute the cfit for each one
for (compound_name in PFASofinterest) {
  # subset the data frame for the current compound
  compound_data <- subset(filter(Baleen_PFAS_data_comb, Sample_type == "baleen"), 
                          Compound == compound_name)
  
  # loop over the ID codes and compute the cfit for each one
  for (id in unique(compound_data$ID_code)) {
    # check if the individual had a baleen sample
    if (id %in% unique(compound_data$ID_code)) {
      # subset the data frame for the current ID code
      id_data <- subset(compound_data, ID_code == id)
      
      # compute the cfit for the current ID code, with a try-catch block to handle errors
      tryCatch(
        {
          mycenfit <- cfit_mod(id_data$Conc_NoCorr + 0.0000001, id_data$Censor_TF)
          # add the cfit result to the list with the compound name and ID code as the name
          mycenfit_list_ID[[paste(compound_name, id, sep="_")]] <- mycenfit
        },
        error = function(e) {
          cat(paste("Error:", e$message, "Skipping individual", id, "for compound", compound_name, "\n"))
        }
      )
    } else {
      cat(paste("Individual", id, "did not have a baleen sample for compound", compound_name, "\n"))
    }
  }
}





# create an empty list to store the mean and sd values
mean_sd_list <- list()

# create an empty data frame to store the results
mycenfit_df <- data.frame(Compound = character(),
                          ID_code = character(),
                          point_est = numeric(),
                          ci_lwr = numeric(),
                          ci_upr = numeric(),
                          stringsAsFactors = FALSE)

# loop over the compound and ID combinations and extract the point estimate, lower and upper bounds
for (compound_id in names(mycenfit_list_ID)) {
  # extract the values from the cfit output
  values <- unlist(lapply(mycenfit_list_ID[[compound_id]], function(x) c(x$point.est, x$ci.lwr, x$ci.upr)))
  
  # add the values to the data frame with the compound and ID as separate columns
  mycenfit_df[nrow(mycenfit_df) + 1,] <- c(strsplit(compound_id, "_")[[1]], values)
}

# convert the columns to their appropriate data types
mycenfit_df$point_est <- as.numeric(mycenfit_df$point_est)
mycenfit_df$ci_lwr <- as.numeric(mycenfit_df$ci_lwr)
mycenfit_df$ci_upr <- as.numeric(mycenfit_df$ci_upr)




for (compound_id in names(mycenfit_list_ID)) {
  # print the element to see what it looks like
  print(mycenfit_list_ID[[compound_id]])
  
  # check if the list element contains the expected structure
  if (!is.na(mycenfit_list_ID[[compound_id]]$point.est)) {
    # extract the values from the cfit output
    values <- unlist(lapply(mycenfit_list_ID[[compound_id]], function(x) c(x$point.est, x$ci.lwr, x$ci.upr)))
    # add the values to the list with the compound and ID as the name
    values_list[[compound_id]] <- values
  } else {
    cat(paste("No valid cfit output for", compound_id, "\n"))
  }
}





# create an empty list to store the results
mycenfit_df_list <- list()

# loop over the compound and ID combinations and extract the three values
for (compound_id in names(mycenfit_list_ID)) {
  # extract the three values from the cfit output
  values <- unlist(lapply(mycenfit_list_ID[[compound_id]], function(x) x[c("point.est", "ci.lwr", "ci.upr")]))
  # add the values to the list with the compound and ID as the name
  mycenfit_df_list[[compound_id]] <- as.data.frame(t(values))
}

# combine all the dataframes in the list into one dataframe
mycenfit_df <- bind_rows(mycenfit_df_list, .id = "compound_id")







# convert the list to a data frame in long tidy format
mean_sd_df <- as.data.frame(do.call(rbind, mean_sd_list), stringsAsFactors = FALSE)
# add column headers
colnames(mean_sd_df) <- c("mean", "sd")
# add Compound and ID_code columns by splitting the names column
mean_sd_df$Compound <- sapply(strsplit(row.names(mean_sd_df), "_"), `[`, 1)
mean_sd_df$ID_code <- sapply(strsplit(row.names(mean_sd_df), "_"), `[`, 2)
# reorder the columns
mean_sd_df <- mean_sd_df[, c("Compound", "ID_code", "mean", "sd")]







# loop over the names in the mycenfit_list_ID and extract the quantiles for each one
for (name in names(mycenfit_list_ID)) {
  # extract the ID_code and Compound name from the name string
  ID_code <- strsplit(name, "_")[[1]][1]
  Compound <- strsplit(name, "_")[[1]][2]
  
  # extract the values from the mycenfit object if it is not empty
  if (!is.null(mycenfit_list_ID[[name]])) {
    # extract the values and force them to be numeric
    values <- as.numeric(mycenfit_list_ID[[name]]$Q50)
    
    # print the values, replacing non-numeric values with "NA"
    cat(paste(ID_code, Compound, paste(ifelse(is.na(values) | !is.numeric(values), "NA", values), collapse=", "), "\n"))
  }
}









