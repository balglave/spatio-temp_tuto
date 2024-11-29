#-----------------------------------------------------
## Build matrix/tensors for the multi-species datasets
#-----------------------------------------------------

Ztensor_full <- array(NA,dim = c(nrow(loc_x),max(time.step_df$t),length(spp_analysis_vec)))
for(i in 1:length(spp_analysis_vec)){
  spp_i <- spp_analysis_vec[i]
  print(paste0("Species: ", spp_i," | Zone: ",area_i))
  
  ## Filter on area and species
  S_x_df_full_2 <- S_x_df_full %>%
    filter(Area == area_i) %>% 
    filter(Species %in% spp_i)
  
  ## Build the spatio-temporal matrix
  S_x_df_full_3 <- S_x_df_full_2 %>%
    mutate(S_x = log(S_x)) %>% 
    dplyr::select(cell,S_x,Year_Month) %>%
    distinct(cell,S_x,Year_Month,keep_all=T) %>%
    arrange(cell) %>% 
    dplyr::select(-keep_all) %>% 
    pivot_wider(values_from = S_x, names_from = cell) %>% 
    dplyr::select(-Year_Month)
  EOF_mat <- as.matrix(S_x_df_full_3)
  spat_mean <- apply(EOF_mat, 2, mean)
  
  ## Center raw data and compute anomaly matrix
  nT <- nrow(EOF_mat)
  Zspat_detrend <- EOF_mat - outer(rep(1, nT), spat_mean)
  Zt <- Zspat_detrend
  Zs <- t(Zt)
  
  ## Make full EOFs matrices
  if(spp_i == spp_analysis_vec[1] & area_i == area_vec[1]){
    Zs_full <- Zs2_full <- Zs
  }else{
    Zs_full <- cbind(Zs_full,Zs)
    Zs2_full <- rbind(Zs2_full,Zs)
  }
  
  ## Fill the tensor
  Ztensor_full[,,i] <- Zs
  
}
