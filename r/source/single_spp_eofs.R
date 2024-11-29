#----------------------
## Single variable EOFs
#----------------------

## Filter on area and species
S_x_df_full_2 <- S_x_df_full %>%
  filter(Area == area_i) %>% 
  filter(Species %in% spp_i)

## Build the spatio-temporal matrix
S_x_df_full_3 <- S_x_df_full_2 %>%
  # mutate(S_x = log(S_x)) %>% # pass to log or not 
  dplyr::select(cell,S_x,Year_Month) %>%
  distinct(cell,S_x,Year_Month,keep_all=T) %>%
  arrange(cell) %>% 
  dplyr::select(-keep_all) %>% 
  pivot_wider(values_from = S_x, names_from = cell) %>% 
  dplyr::select(-Year_Month)

EOF_mat <- as.matrix(S_x_df_full_3)

spat_mean <- apply(EOF_mat, 2, mean) # spatial average pattern
temp_mean <- apply(EOF_mat, 1, mean) # temporal average pattern

## Center raw data
nT <- nrow(EOF_mat)
nS <- ncol(EOF_mat)
Zspat_detrend <- EOF_mat - outer(rep(1, nT), spat_mean) # center on spatial average pattern
temp_mean <- apply(Zspat_detrend, 1, mean) # temporal average pattern
spat_mean <- apply(Zspat_detrend, 2, mean) # spatial average pattern
Zspat_detrend <- t(t(Zspat_detrend) - outer(rep(1, nS), temp_mean)) # center on spatial average pattern

## Compute variance ...
s2_t <- apply(t(EOF_mat),2,"var",na.rm=T) # ... relatively to the space (results in a time series of variance)
s2_s <- apply(t(EOF_mat),1,"var",na.rm=T) # ... relatively to the time (results in a map of variance)

## Reduce the matrices or not
Zt <- Zspat_detrend # %*% diag(1/sqrt(s2_s))
Zs <- t(Zt) # %*% diag(1/sqrt(s2_t))

## Perform svd
E <- svd((Zs)) # svd

sum(apply(E$u[,1:2] %*% diag(E$d[1:2]) %*% t(E$v[,1:2]),2,var)) /
  sum(apply((Zs),2,var))

sum((E$d^2 / sum(E$d^2))[1:2])

n_EOF <- 2 # number of dimension to retain in the analysis

## Spatial patterns
SpatPat <- E$u
colnames(SpatPat) <- paste0("Dim", 1:nrow(EOF_mat)) # label columns
EOFs <- cbind(loc_x[,c("long","lati","cell")], SpatPat)
EOFs_2 <- EOFs %>%
  select_at(c("long","lati","cell",paste0("Dim",1:n_EOF))) %>%
  pivot_longer(cols = starts_with("Dim")) %>%
  dplyr::rename(EOF = name) %>%
  filter(EOF %in% c(paste0("Dim",1:n_EOF)))

EOFs_plot = ggplot()+
  geom_point(data=EOFs_2,(aes(x=long,y=lati, col = -value)))+
  scale_color_continuous_divergingx(palette = "Spectral", mid = 0,rev = T)+
  theme_bw()+
  xlab("") + ylab("")+
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90))+
  facet_wrap(.~EOF)

## Temporal patterns
EOFset1 <- E$v[1:(nT),]
EOFset1 <- EOFset1 %>% data.frame
colnames(EOFset1) <- paste0("Dim",1:ncol(EOFset1))
EOFset2 <- EOFset1[,1:n_EOF]
EOFset1_2 <- EOFset2 %>%
  mutate(t = 1:nrow(EOFset2)) %>%
  pivot_longer(cols = starts_with("Dim")) %>%
  dplyr::rename(PC = name) %>%
  filter(PC %in% c(paste0("Dim",1:n_EOF))) %>%
  inner_join(time.step_df)

EOFt_plot = ggplot(EOFset1_2)+
  geom_hline(yintercept=0,linetype="dashed", color = "black", size = 0.2)+
  geom_line(size=1,aes(x=Year_Month,y=value,col=PC,group=PC))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),aspect.ratio = c(1/4))+
  xlab("")+ylab("")+
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))+
  coord_cartesian(clip = "off")+
  facet_wrap(.~PC,ncol = 1)+
  geom_vline(xintercept = EOFset1_2$Year_Month[which(str_detect(EOFset1_2$Year_Month,"01-01"))],
             linetype="dashed",color = "skyblue")
