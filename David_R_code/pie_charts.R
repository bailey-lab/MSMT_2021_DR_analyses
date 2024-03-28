######################
#kagera 561H district map
######################

#load required libraries for base maps
library(tidyverse)
library(rnaturalearth)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(scatterpie)

rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines',
                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "medium", type = 'lakes',
                       category = 'physical', returnclass = "sf")
oceans10 <- ne_download(scale = "medium", type = "coastline",
                        category = 'physical', returnclass = "sf")
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")
admin2 <- st_read("tza_admbnda_adm2_20181019.shp")
admin10 <- st_read("tza_admbnda_adm1_20181019.shp")
#read facility data
kagera_aatable_bin_561_health_facility <- read.csv("R561H_facility_final.csv")
#make percent
kagera_aatable_bin_561_health_facility$R561H <- kagera_aatable_bin_561_health_facility$R561H*100
kagera_aatable_bin_561_health_facility$WT <- 100-kagera_aatable_bin_561_health_facility$R561H
ggplot()+
  geom_sf(data=lakes10, color="grey95", fill ="aliceblue", size= 0.8) +
  geom_sf(data=admin10, color="grey40", size= 1, alpha = 0.1, fill = "white") +
  geom_sf(data=admin2, color="grey40", size= 1, alpha = 0.1, fill = "grey60") +
  geom_scatterpie(aes(x=lon, y=lat), data=kagera_aatable_bin_561_health_facility,
                  cols = c("R561H", "WT"),
                  pie_scale = 4, color = "grey20", alpha=.8, legend_name = "K13") +
  scale_fill_manual(values = c("#AC2376B3", "white")) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.85, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  annotate("text", x = 30.62, y = -3.12, label = "Burundi",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.6, y = -1.7, label = "Rwanda",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.9, y = -0.9, label = "Uganda",
           color="grey60", size=4 , fontface="italic") +
  coord_sf(xlim = c(30.3, 32), ylim = c(-0.89,-3.4), expand = TRUE) +
  theme(plot.background = element_rect(fill = "white"))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12)) +
  theme_void()

#save plot
ggsave("/home/dgiesbre/Documents/bailey_pubs/MSMT_DR_2021/kagera_R561H.png", dpi=600, width=6, height=6)
ggsave("/home/dgiesbre/Documents/bailey_pubs/MSMT_DR_2021/kagera_R561H.svg", dpi=600, width=6, height=6)