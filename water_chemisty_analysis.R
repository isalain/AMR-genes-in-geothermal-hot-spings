# Create the a plot of altitude and water temperature differences
alt_temper<-ggscatterhist(
  site_data, x = "Temperature", y = "Elevation",
  color = "Location", size = 3, alpha = 0.6,
  palette = c("#009E73", "#E7B800"),
  margin.plot = "histogram",
  ggtheme = theme_bw())
# Comparison of water chemistry in Bugarama and Gisenyi GHS 
temp_boxplot<- ggviolin(site_data, "Location", "Temperature",
                        fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")
#pH
pH_boxplot<- ggviolin(site_data, "Location", "pH",
                        fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#TDS
TDS_boxplot<- ggviolin(site_data, "Location", "TDS",
                      fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#EC
EC_boxplot<- ggviolin(site_data, "Location", "EC",
                       fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#ORP
ORP_boxplot<- ggviolin(site_data, "Location", "ORP",
                      fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#TP
TP_boxplot<- ggviolin(site_data, "Location", "TP",
                       fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#TN
TN_boxplot<- ggviolin(site_data, "Location", "TN",
                      fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#NO3N
NO3N_boxplot<- ggviolin(site_data, "Location", "NO3N",
                      fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")
#NH3N
NH3N_boxplot<- ggviolin(site_data, "Location", "NH3N",
                        fill = "Location", palette = c("#009E73", "#E7B800"))+ 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3)+
  stat_kruskal_test()+xlab("")+
  theme(legend.position = "none")





