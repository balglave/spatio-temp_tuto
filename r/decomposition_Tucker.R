library(tidyverse)
library(rTensor)
library(colorspace)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)

## cell = cellule
## s_x = prediction (champ latent)
## area = zone géograhique (mer celtique, manche et gascogne)

#-----------------------------------------------------
## Démarche
#-----------------------------------------------------

# 1. Faire EOF sur tenseur
# 2. Faire clustering des espèces pour regrouper celles qui sont similaires
# 3. Faire une seule carte et un seul loading par groupe d'espèces

load("species_datarmor.RData")

S_x_df_full$Species <- as.factor(S_x_df_full$Species)

area_vec <- c("bob_cs")

S_x_df_full_2 <- S_x_df_full %>%
  filter(Area == area_vec)

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

#-----------------------------------------------------
## Build tensors for the multi-species datasets
#-----------------------------------------------------
S_x_df_full_2$Species <- droplevels(S_x_df_full_2$Species)
spp_analysis_vec <- levels(S_x_df_full_2$Species)
spp_analysis_vec <- spp_analysis_vec[spp_analysis_vec != "Raja_brachyura"]
print(spp_analysis_vec)

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

save(Ztensor_full,
     file = "data/tensor_19_spp.Rdata")

#-----------------------------------------------------
## Tucker decomposition
#-----------------------------------------------------

load("data/tensor_19_spp.Rdata")

graph_Tucker <- function(YY,liste_tps = NULL, dim=1:2){
  if (!is.null(liste_tps)) YY <- YY[,liste_tps,]
  YY_ts <- rTensor::as.tensor(YY)
  res_tucker <- rTensor::hosvd(YY_ts,ranks=c(max(dim),max(dim),max(dim))) ## décomposition de Tucker
  
  ## Liste de matrices factorielles : 
  ## - res_tucker$U[[1]] = facteurs spatiaux
  ## - res_tucker$U[[2]] = facteurs temporels
  ## - res_tucker$U[[3]] = facteurs liés aux espèces
  
  gr1 <- data.frame(loc_x,spatial_factor = res_tucker$U[[1]][,dim[1]]) %>% 
    ggplot()+  geom_point(aes(x=long,y=lati,col=spatial_factor))+
    ggtitle(paste0("Dim",dim[1]))+
    scale_fill_continuous_divergingx(palette = "Spectral", mid = 0,rev = T)+ 
    scale_color_continuous_divergingx(palette = "Spectral", mid = 0,rev = T)+coord_fixed(ratio = 1) ## dimensions 1
  gr2 <- data.frame(loc_x,spatial_factor = res_tucker$U[[1]][,dim[2]]) %>% 
    ggplot()+  geom_point(aes(x=long,y=lati,col=spatial_factor))+
    ggtitle(paste0("Dim",dim[2]))+
    scale_fill_continuous_divergingx(palette = "Spectral", mid = 0,rev = T)+ 
    scale_color_continuous_divergingx(palette = "Spectral", mid = 0,rev = T)+coord_fixed(ratio = 1) ## dimension 2
  
  gr1t <- data.frame(t=1:nrow(res_tucker$U[[2]]),loadings = res_tucker$U[[2]][,dim[1]]) %>% 
    ggplot(aes(x=t,y=loadings))+   geom_line() + ggtitle(paste0("Dim",dim[1]," - temp")) ## dimensions 1
  gr2t <- data.frame(t=1:nrow(res_tucker$U[[2]]),loadings = res_tucker$U[[2]][,dim[2]]) %>% 
    ggplot(aes(x=t,y=loadings))+  geom_line() + ggtitle(paste0("Dim",dim[2]," - temp")) ## dimensions 2
  
  df <- data.frame(tps = 1:nrow(res_tucker$U[[2]]), dim1 = res_tucker$U[[2]][,dim[1]], dim2 = res_tucker$U[[2]][,dim[2]])
  gr3 <- ggplot(df, aes(x = dim1, y = dim2, color = tps)) +
    geom_path() +  scale_color_gradient2(low = "blue",mid="yellow",midpoint =mean(res_tucker$U[[2]][,dim]),  high = "red") + 
    theme_minimal()+coord_fixed(ratio = 1) + ggtitle("Coeff du temps")+xlab(paste0("Dim",dim[1]))+ylab(paste0("Dim",dim[2]))
  
  group <- res_tucker$U[[3]][,dim]
  gr4 <- data.frame(group,esp=rep(spp_analysis_vec,nrow(group)/19),an=rep(1:(nrow(group)/19),each=19)) %>%
    ggplot()+coord_fixed(ratio = 1) +  aes(x=group[,1],y=group[,2],col=esp,group=esp)+
    geom_point()+xlab(paste0("Dim",dim[1]))+ylab(paste0("Dim",dim[2]))+ggtitle("Représentation espèce")
  
  (gr1+gr2)/(gr3+gr4)
}

ZZ <- Ztensor_full

## Standardisation (centering and reduction)
for (i in 1:dim(ZZ)[3]) {
  #  Ztensor_full[,,i] <- Ztensor_full[,,i]-mean(as.vector(Ztensor_full[,,i]))
  #  Ztensor_full[,,i] <- Ztensor_full[,,i]/sqrt(var(as.vector(Ztensor_full[,,i])))
  ZZ[,,i] <- sweep(ZZ[,,i],2,colMeans(ZZ[,,i]),FUN="-")
  ZZ[,,i] <- sweep(ZZ[,,i],2,apply(ZZ[,,i],2,sd),FUN="/")
}

## Plot Tucker svd outputs
tucker_plot <- graph_Tucker(ZZ,liste_tps = 1:180, dim=c(1,2))
plot(tucker_plot)


#-----------------------------------------------------
## Clustering espèces
#-----------------------------------------------------

dim=1:2
YY_ts <- rTensor::as.tensor(ZZ)
res_tucker <- rTensor::hosvd(YY_ts,ranks=c(max(dim),max(dim),max(dim))) ## décomposition de Tucker
species_factors <- res_tucker$U[[3]]

res_tucker$est

dist_matrix <- dist(species_factors, method = "euclidean")

hcpc <- hclust(dist_matrix, method = "ward.D2")
plot(hcpc, main = "Dendrogramme des espèces", xlab = "Espèces", ylab = "Distance")
clusters <- cutree(hcpc, k = 4)
table(clusters)

# Ajout des clusters au dataframe gr4
group <- res_tucker$U[[3]][,dim]
gr4 <- data.frame(group, esp = rep(spp_analysis_vec, nrow(group) / 19), 
                  an = rep(1:(nrow(group) / 19), each = 19), cluster = as.factor(clusters),
                  species_name = spp_analysis_vec)

ggplot(gr4) + 
  coord_fixed(ratio = 1) +
  aes(x = group[, 1], y = group[, 2], col = cluster, group = cluster) + 
  geom_point() + 
  geom_text_repel(aes(label = species_name), size = 3) +
  xlab(paste0("Dim", dim[1])) + 
  ylab(paste0("Dim", dim[2])) + 
  ggtitle("Représentation des espèces avec clusters")
save(res_tucker,
     file = "data/res_tucker.Rdata")

#-----------------------------------------------------
## Pattern spatial par cluster d'espèces
#-----------------------------------------------------

## Reconstruction of partial tensor for cluster 1

for (i in 1:4){
  index_cluster <- which(clusters == i)
  U3_cluster <- res_tucker$U[[3]][index_cluster,]
  ZZ_clust <- ttl(res_tucker$Z, list(res_tucker$U[[1]], res_tucker$U[[2]], U3_cluster), ms = c(1,2,3))
}


index_cluster1 <- which(clusters == 1)
U3_cluster1 <- res_tucker$U[[3]][index_cluster1,]
ZZ_clust1 <- ttl(res_tucker$Z, list(res_tucker$U[[1]], res_tucker$U[[2]], U3_cluster1), ms = c(1,2,3))
dim(ZZ_clust1)