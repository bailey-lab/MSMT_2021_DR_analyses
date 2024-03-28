################################################
#prevalence chloropleths
################################################
#load required libraries for base maps
library(tidyverse)
library(rnaturalearth) 
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
#2021 prevalence data clinical variants

rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', 
                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "medium", type = 'lakes', 
                       category = 'physical', returnclass = "sf") #useful to have these two for points of reference
oceans10 <- ne_download(scale = "medium", type = "coastline",
                        category = 'physical', returnclass = "sf")
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")
admin10 <- st_read("tza_admbnda_adm1_20181019.shp") #can get a pretty good map just with this
                        category = 'physical', returnclass = "sf")
##########################################
#561H national plot
##########################################
clinical_variants <- read.csv("231106_dhps-dhfr_regional.csv")

tanz_regions <- dplyr::right_join(admin10,clinical_variants, by = "ADM1_EN" )
#make percent
tanz_regions$k13_R561H<- tanz_regions$k13_Arg561His * 100
# Create breaks and discretize values
br <- c(-0.1,0,2,4,6,8)
#br <- c(-0.1,0,1,2,4,6,8)
tanz_regions$k13_R561H<- cut(tanz_regions$k13_R561H,
                             breaks = br,
                             dig.lab = 5)
summary(tanz_regions$k13_R561H)
# Palettes
pal <- hcl.colors(9, "SunsetDark", rev = TRUE, alpha = 0.7)
pal2 <- hcl.colors(7, "TealRose", alpha = 0.7)
palzero <- c("#DADEDF","#FFD99FB3" ,"#FFAC7CB3", "#FB796CB3", "#E74075B3","#B52578B3","#B52578B3")
#custom labels
labs <- c("0%","0-2%", "2-4%", "4-6%", "6-8%")
#plot figure 2A
ggplot()+
  geom_sf(data=sov110, color='grey30', size=0.8, alpha = 0.2, linewidth = 1,
          fill = ifelse(sov110$ADMIN == "United Republic of Tanzania", 'white', 'seagreen1')) +
  geom_sf(data = tanz_regions,
          aes(fill = k13_R561H), color = NA) +
  scale_fill_manual(values = palzero, drop = FALSE, na.value = "grey96",
                    labels = labs, na.translate = FALSE, name= "k13_R561H\nPrevalence",
                    guide = guide_legend(direction = "vertical",
                                         nrow = 11, label.position = "right", 
                                         name = "NULL")) +

  geom_sf(data=lakes10, color="grey40", fill ="lightblue", size= 0.8) +
  annotate("text", x = 33.8, y = -11, label = "Malawi",
           color="grey60", size=3 , fontface="italic") +
  annotate("text", x = 32, y = -10, label = "Zambia",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 38, y = -12.0, label = "Mozambique",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 38, y = -2.5, label = "Kenya",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30, y = -2.0, label = "Rwanda",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 29.5, y = -7.3, label = "DRC",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30, y = -3.2, label = "Burundi",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.9, y = -0.7, label = "Uganda",
           color="grey60", size=4 , fontface="italic") +
  geom_sf(data=admin10, color="grey40", size= 0.6, alpha = 0.1) +
  coord_sf(xlim = c(29.5, 40.2), ylim = c(-1,-12), expand = TRUE) +
  theme_void()+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))
    
      
#save plot
ggsave("fig2A_k13_R561H.svg", dpi=600, width=6, height=5)
ggsave("fig2A_k13_R561H.png", dpi=600, width=6, height=5)

###########################
#581
###########################
clinical_variants <- read.csv("231106_dhps-dhfr_regional.csv")

tanz_regions <- dplyr::right_join(admin10,clinical_variants, by = "ADM1_EN" )
#make percent
tanz_regions$dhps_A581G <- tanz_regions$dhps_Ala581Gly * 100
# Create breaks and discretize values
br <- c(-0.1,10,20,30,40,50)
#br <- c(-0.1,0,1,2,4,6,8)
tanz_regions$dhps_A581G <- cut(tanz_regions$dhps_A581G,
                              breaks = br,
                              dig.lab = 5)
summary(tanz_regions$dhps_A581G)
# Palettes
pal <- hcl.colors(6, "SunsetDark", rev = TRUE, alpha = 0.7)
pal2 <- hcl.colors(7, "TealRose", alpha = 0.7)
palzero <- c("#DADEDF","#FFD99FB3" ,"#FFC38DB3", "#FFAC7CB3", "#FF9471B3", "#FB796CB3" ,"#F65A6DB3", "#E74075B3", "#CF3179B3" ,"#B52578B3","#7D1D67B3")
#custom labels
#labs <- c(-0.1,0,1,2,4,6,8,10)
labs <- c("0-9%","10-20%", "20-30%", "30-40%", "40-50%","50-60%")
#plot map fig. 2D
ggplot()+
  geom_sf(data=sov110, color='grey30', size=0.8, alpha = 0.2, linewidth = 1,
          fill = ifelse(sov110$ADMIN == "United Republic of Tanzania", 'white', 'seagreen1')) +
  geom_sf(data = tanz_regions,
          aes(fill = dhps_A581G), color = NA) +
  scale_fill_manual(values = pal, drop = FALSE, na.value = "grey96",
                    labels = labs, na.translate = FALSE, name= "dhps_A581G\nPrevalence",
                    guide = guide_legend(direction = "vertical",
                                         nrow = 11, label.position = "right", 
                                         name = "NULL")) +
  geom_sf(data=lakes10, color="grey40", fill ="lightblue", size= 0.8) +
  annotate("text", x = 33.8, y = -11, label = "Malawi",
           color="grey60", size=3 , fontface="italic") +
  annotate("text", x = 32, y = -10, label = "Zambia",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 38, y = -12.0, label = "Mozambique",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 38, y = -2.5, label = "Kenya",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30, y = -2.0, label = "Rwanda",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 29.5, y = -7.3, label = "DRC",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30, y = -3.2, label = "Burundi",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.9, y = -0.7, label = "Uganda",
           color="grey60", size=4 , fontface="italic") +
  geom_sf(data=admin10, color="grey40", size= 0.6, alpha = 0.1) +
  coord_sf(xlim = c(29.5, 40.2), ylim = c(-1,-12), expand = TRUE) +
  theme_void()+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))
#save plot
ggsave("fig2d_dhps_581.svg", dpi=600, width=6, height=5)
ggsave("fig2d_dhps_581.png", dpi=600, width=6, height=5)




##########################################
#same plot for 164L
##########################################
admin10 <- st_read("tza_admbnda_adm1_20181019.shp") #can get a pretty good map just with this
category = 'physical', returnclass = "sf")

tanz_regions <- dplyr::right_join(admin10, clinical_variants, by = "ADM1_EN" )
tanz_regions$dhfr_Ile164Leu
#make percent
tanz_regions$dhfr_I164L <- tanz_regions$dhfr_Ile164Leu * 100
tanz_regions$dhfr_I164L
# Create breaks and discretize values
br <- c(-0.1,0,2,4,6,8,10,12,14,16)
#br <- c(-0.1,0,1,2,4,6,8)
tanz_regions$dhfr_I164L <- cut(tanz_regions$dhfr_I164L,
                               breaks = br,
                               dig.lab = 5)
summary(tanz_regions$dhfr_I164L)
# Palettes
pal <- hcl.colors(9, "SunsetDark", rev = TRUE, alpha = 0.7)
pal2 <- hcl.colors(7, "TealRose", alpha = 0.7)
palzero <- c("#DADEDF","#FFD99FB3" ,"#FFC38DB3", "#FFAC7CB3", "#FF9471B3", "#FB796CB3" ,"#F65A6DB3", "#E74075B3", "#CF3179B3" ,"#B52578B3","#B52578B3")
#custom labels
labs <- c("0%","0-2%", "2-4%", "4-6%", "6-8%","8-10%", "10-12%", "12-14%", "14-16%")
#plot map fig2C_dhfr
ggplot()+
  geom_sf(data=sov110, color='grey30', size=0.8, alpha = 0.2, linewidth = 1,
          fill = ifelse(sov110$ADMIN == "United Republic of Tanzania", 'white', 'seagreen1')) +
  geom_sf(data = tanz_regions,
          aes(fill = dhfr_I164L), color = NA) +
  scale_fill_manual(values = palzero, drop = FALSE, na.value = "grey96",
                    labels = labs, na.translate = FALSE, name= "dhfr_I164L\nPrevalence",
                    guide = guide_legend(direction = "vertical",
                                         nrow = 11, label.position = "right", 
                                         name = "NULL")) +
  geom_sf(data=lakes10, color="grey40", fill ="lightblue", size= 0.8) +
  annotate("text", x = 33.8, y = -11, label = "Malawi",
           color="grey60", size=3 , fontface="italic") +
  annotate("text", x = 32, y = -10, label = "Zambia",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 38, y = -12.0, label = "Mozambique",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 38, y = -2.5, label = "Kenya",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30, y = -2.0, label = "Rwanda",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 29.5, y = -7.3, label = "DRC",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30, y = -3.2, label = "Burundi",
           color="grey60", size=4 , fontface="italic") +
  annotate("text", x = 30.9, y = -0.7, label = "Uganda",
           color="grey60", size=4 , fontface="italic") +
  geom_sf(data=admin10, color="grey40", size= 0.6, alpha = 0.1) +
  coord_sf(xlim = c(29.5, 40.2), ylim = c(-1,-12), expand = TRUE) +
  theme_void()+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))
#save plot
ggsave("fig2c_dhfr_I164L.svg", dpi=600, width=6, height=5)
ggsave("fig2c_dhfr_I164L.png", dpi=600, width=6, height=5)

