#---------------------------------------------------------------------
## Workflow for analysing spatio-temporal data of species distribution
#---------------------------------------------------------------------

library(colorspace)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

load("data/S_x_df_full.RData")

spp_analysis_vec <- c("Solea_solea","Dicentrarchus_labrax","Sepia_officinalis") # species under study
area_vec <- c("bob_cs") # area under study

loc_x <- S_x_df_full %>% 
  filter(Species %in% spp_analysis_vec) %>% 
  filter(Area %in% area_vec) %>% 
  group_by(long,lati,cell) %>% 
  slice(1) %>% 
  dplyr::select(long,lati,cell) %>% 
  arrange(cell)

time.step_df <- S_x_df_full %>% 
  filter(Species %in% spp_analysis_vec) %>% 
  filter(Area %in% area_vec) %>% 
  group_by(Month,Year,Year_Month,t) %>% 
  slice(1) %>% 
  dplyr::select(Month,Year,Year_Month,t) %>% 
  arrange(t)

#-----------------------------------
## Perform EOFs for a single species
#-----------------------------------
spp_i <- "Solea_solea"
area_i <- "bob_cs"
source("r/source/single_spp_eofs.R")

#-----------------------------------------------------
## Build matrix/tensors for the multi-species datasets
#-----------------------------------------------------
source("r/source/build_matrix_multi.R")

