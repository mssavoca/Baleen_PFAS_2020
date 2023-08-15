--------------------------------------------------------------------------------
  # Stats for PFAS-baleen paper
  # Matt Savoca
  # Started on: 7/28/23
--------------------------------------------------------------------------------
  
  
# not enough data to run the model
#Still neeed to work on  
  
PFAS_baleen_prey_m1 <- lmer(Total_PFAS ~ Prey + (1|ID_code), 
                            data = TL_PFAS_idv)



