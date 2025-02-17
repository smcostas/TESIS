
library(tidyverse)
library(viridis)
library(patchwork)
library(ggcorrplot)
#df_completo <- read.csv('df_completo.csv')
#df_resumen <- read.csv('tabla_resumen.csv')
#colnames(df_resumen);colnames(df_completo)

getwd()      
#setwd("C:/Users/santi/OneDrive/Escritorio/PC-IMBIV/doctorado/asteraceae - muestreo/segundo capitulo/tablas segundo capitulo")
df_completo <- read.csv('df_completo.csv')
df_completo = df_completo[df_completo$spp != 'Thelesperma_megapotamicum',]
df_completo$tiempo <- df_completo$tiempo/1000000
df_completo$sv <- 2.2/df_completo$tiempo
df_completo$log.sv <- log(df_completo$sv)
df_completo$log.masa <- log(df_completo$masa_mg)
df_completo$area_p = pi*df_completo$l_papus^2

max(df_completo$sv)
min(df_completo$sv)
sd(df_completo$sv)
mean(df_completo$sv)
plot(df_completo$log.sv)
head(df_completo)

p_sv <- ggplot(df_completo, aes(x = sv))
p_sv <- p_sv + geom_histogram(aes(x = sv, y = after_stat(density)),
                                    col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                                    fill = 'white', bins = 30) + scale_x_continuous(breaks = seq(0,3.5, by = 0.5)) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) +
  xlab('Velocidad de Asentamiento (m/s)') + ylab('Densidad') +
  scale_y_continuous(expand = c(0,0))+ 
  theme(
    axis.text.x = element_text(size = 12, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(colour = "grey50"),
    plot.margin = margin(6,0,0,6),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = 'white', color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  )
p_sv




df_completo[df_completo$tipo_papus == 'epappose',]$l_papus = 0
df_completo[df_completo$tipo_papus == 'epappose',]$area_p = 0
df_completo$pl_m = df_completo$masa_mg/(df_completo$area_p+0.0000001)
df_sum <-df_completo[!is.na(df_completo$masa_mg) & !is.na(df_completo$l_papus), ] %>% group_by(spp)  %>% 
  summarise_at(vars(sv, masa_mg, log.sv, log.masa, l_papus, area_p, pl_m), list(mean))
colnames(df_sum) <- c('spp', 'sv', 'masa', 'log.sv', 'log.masa', 'l_papus', 'area_p', 'pl_m')
colnames(df_completo)
empty_col <- c('tipo_papus','Latitud','Longitud','altitud','Ecoregion','habito','muestreo','Tribu', 'cod_gps')
df_sum[,empty_col] <- NA
for (i in 1:76){
  for (c in 1:615){
    if (df_sum$spp[i] == df_completo$spp[c]){
      df_sum$tipo_papus[i] <- df_completo$tipo_papus[c]
      df_sum$Latitud[i] <- df_completo$Latitud[c]
      df_sum$Longitud[i] <- df_completo$Longitud[c]
      df_sum$altitud[i] <- df_completo$altitud[c]
      df_sum$Ecoregion[i] <- df_completo$Ecoregion[c]
      df_sum$habito[i] <- df_completo$habito[c]
      df_sum$muestreo[i] <- df_completo$muestreo[c]
      df_sum$Tribu[i] <- df_completo$Tribu[c]
      df_sum$cod_gps[i] <- df_completo$cod_gps[c]
      break
    }
  }
}
library(stringr)

## grafico resumen muestreo #######
subdf <- data.frame(tribe = rep('', length(levels(df_sum$Tribu))), 
                    value = rep(0, length(levels(df_sum$Tribu)))
)

for (i in 1:length(levels(df_sum$Tribu))){
  subdf$tribe[i] <- levels(df_sum$Tribu)[i]
  subdf$value[i] <- length(which(df_sum$Tribu == subdf$tribe[i]))
}
l <- reorder(subdf$tribe, subdf$value , decreasing = T)

subdf$tribe <- factor(subdf$tribe, levels = levels(l))

my_theme <-   theme(
  legend.text = element_text(size = 8),
  legend.background = element_rect(fill = "white"),
  plot.title = element_text(vjust = -1, hjust = 0.5),
  plot.subtitle = element_text(vjust = -1, hjust = 0.5),
  axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
  axis.text.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  axis.line = element_line(colour = "grey50"),
  plot.margin = margin(6,0,0,6),
  panel.grid = element_line(color = "#b4aea9"),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(linetype = "dashed"),
  panel.background = element_rect(fill = 'white', color = "white"),
  plot.background = element_rect(fill = "white", color = "white")
)
bar.plot2 <- ggplot(subdf, aes(x = '', y = value, fill = tribe))
bar.plot2 + geom_bar(stat = 'identity', width=1) + coord_polar('y', start = 0) + 
  scale_fill_viridis_d('Tribes', direction = -1, option = 'A') + theme_void()

bar.plot2 + geom_bar(stat = 'identity') + coord_polar('y', start = 0) + 
  scale_fill_viridis_d('Tribes', option = 'A', direction = -1) + labs(title = 'Tribes Sampled') +
  theme_void() + theme(legend.text = element_text(size = 8),
                       legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                       legend.position = 'bottom',
                       legend.direction = 'horizontal',
                       plot.title = element_text(vjust = -1, hjust = 0.5),
                       panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                       plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
  )





### Ecoregiones ####
### Barplot - distribucion por tribus ####

eco_count <- count(df_sum, Ecoregion)
df_sum$Ecoregion <- factor(df_sum$Ecoregion, 
                       levels = eco_count$Ecoregion[order(eco_count$n, decreasing = T)])

df_sum$Ecoregion = factor(df_sum$Ecoregion, labels = c('Puna andina central', 
                                                       'Monte', 'Chaco seco', 
                                                       'Yungas', 'Estepa patagónica',
                                                       'Bosques pat/estepa subandina'))
bar.plot <- ggplot(df_sum, aes(x = Ecoregion, fill = Tribu))

#png('ecoreg_trib.png', width = 18, height = 18, units = 'cm', res = 300)
#pdf('ecoreg_trib.pdf')
bar.plot + geom_bar() + scale_y_continuous(expand = c(0,0), 
                                           breaks = seq(0,30, by = 10),
                                           limits = c(0,30)) + 
  scale_fill_viridis_d('Tribus', option = 'A', direction = -1) +xlab('Ecoregión') + ylab('n de especies') + 
  my_theme + theme(axis.text.x = element_text(angle = 20, hjust = 0.65, vjust = 0.8))
dev.off()


bar.plot2 <- ggplot(df_sum, aes(x = Tribu, fill = Tribu))
#png('trib.png', width = 18, height = 18, units = 'cm', res = 300)
#pdf('trib.pdf')
bar.plot2 + geom_bar(show.legend = F) + scale_y_continuous(expand = c(0,0), 
                                  breaks = seq(0,15, by = 3),
                                  limits = c(0,15)) + 
  scale_fill_viridis_d('Tribus', option = 'A', direction = -1) +xlab('Tribus') + ylab('n de especies') + 
  my_theme + theme(axis.text.x = element_text(angle = 50, hjust = 0.65, vjust = 0.7))
dev.off()
## continua #######
## para el arbol
sample_species_list = df_sum[,'spp']
sample_species_list$spp = gsub('_', ' ',sample_species_list$spp)


sample_species_list <- sample_species_list %>%
  mutate(genus = word(spp, 1))

#write.csv(sample_species_list, 'sample_species_list_asteraceae.csv')

p <- ggplot(df_sum, aes(x = sv))
p <- p + geom_histogram(aes(x = sv, y = ..density..),
                        col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                        fill = 'white', bins = 30) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) + scale_x_continuous(breaks = seq(0,3.5, by = 0.5)) +
  xlab('Velocidad de Asentamiento (m/s)') + ylab('Densidad') +
  scale_y_continuous(expand = c(0,0))+ 
  theme(
    axis.text.x = element_text(size = 12, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(colour = "grey50"),
    plot.margin = margin(6,0,0,6),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = 'white', color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  )

p


head(df_sum)
df_sum[(df_sum$habito == 'Subarbusto'),'habito'] <- 'subarbusto'

#lbls <- c('uk', 'Puna', 'Chaco', 'Monte', 'Bosque pat', 'Estepa pat', 'Yungas')
#lvls <- c('Yungas', 'Bosque pat', 'Estepa pat', 'Puna', 'Monte', 'Chaco', 'uk')
df_sum$Ecoregion <- as.factor(df_sum$Ecoregion)
df_sum$Tribu = as.factor(df_sum$Tribu)

length(levels(df_sum$Tribu))
levels(df_sum$Ecoregion)
summary(df_sum$Ecoregion)
#levels(df_sum$Ecoregion) <- lbls
#df_sum$Ecoregion <- factor(df_sum$Ecoregion, levels = lvls)
summary(df_sum$Ecoregion)
p1 <- ggplot(data = df_sum, aes(x = sd_sv, y = sv, color = Ecoregion))
p1 + geom_point() + geom_smooth(method = 'lm') + ggtitle('SV ~ sd') + 
  scale_color_viridis_d() + scale_color_viridis_d(option = 'A', end = 0.90) +
  theme_bw()



p2 <- ggplot(data = df_sum, aes (y = sv, x = Ecoregion, color = Ecoregion))
p2 + geom_boxplot(show.legend = F) + scale_color_viridis_d(option = 'A', end = 0.90) + ggtitle('SV ~ Ecoregion') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, vjust = 0.8, hjust = 0.75))



p3 = ggplot(df_sum, aes(x = tipo_papus, fill = tipo_papus))

p3 + geom_histogram(stat = 'count', show.legend = F) + facet_wrap(~Ecoregion) + xlab('tipo de papus') +
  theme_bw() +  scale_fill_viridis_d(option = 'A', end = 0.90) + theme(axis.text.x = element_text(angle = 25, vjust = 0.8))

df_sum$habito <- factor(as.factor(df_sum$habito), levels = c('subarbusto', 'cojin', 'arbusto', 'hierba'))
p4 <- ggplot(data = df_sum, aes (y = sv, x = habito, color = habito))
p4 + geom_boxplot() + scale_color_viridis_d(option = 'A', end = 0.9) + ggtitle('SV ~ Habito') + theme_bw()


p5 <- ggplot(data = df_sum, aes (x = tipo_papus, fill = tipo_papus))
p5  + geom_histogram(stat = 'count', show.legend = F) + facet_wrap(~habito) + xlab('tipo de papus') +
  theme_bw() +  scale_fill_viridis_d(option = 'A', end = 0.90) +  theme(axis.text.x = element_text(angle = 25, vjust = 0.8))



df_sum$Tribu <- as.factor(df_sum$Tribu)
df_sum$Tribu <- reorder(df_sum$Tribu, df_sum$Tribu, FUN=length, decreasing = T)
p4 <- ggplot(data = df_sum, aes (y = sv, x = Tribu, color = Tribu))
p4 <- p4 + geom_boxplot() + scale_x_discrete(labels = NULL)  + scale_color_viridis_d(option = 'A') + ggtitle('SV ~ Tribu') + theme_bw()

p4
asteroideae <- c('Heliantheae', 'Eupatorieae', 'Tageteae', 'Millerieae', 'Helenieae',
                 'Inuleae', 'Gnaphalieae', 'Astereae', 'Anthemideae', 'Senecioneae',
                 'Calenduleae', 'Bahieae', 'Chaenactidae','Coreopsideae', 'Madieae', 'Perityleae', 'Pertyeae')
cichorioideae <- c('Arctotideae', 'Cichorieae', 'Vernonieae', 'Liabeae','Eremothamneae',
                   'Moquineae', 'Platycarpheae','Athroismeae')

gochnatioideae <- c('Gochnatieae')

mutisioideae <- c('Mutisieae', 'Nassauvieae', 'Onoserideae')

stifftioideae <- c('Hyalideae', 'Stifftieae')

Barnadesioideae <- c('Barnadesieae')

df_sum$subfam <- NA

for (i in 1:76){
  if (df_sum$Tribu[i] %in% asteroideae) {
   df_sum$subfam[i] <- 'Asteroideae'
   }
  else if (df_sum$Tribu[i] %in% cichorioideae) {
    df_sum$subfam[i] <- 'Cichorioideae'
    }
  else if (df_sum$Tribu[i] %in% gochnatioideae) {
    df_sum$subfam[i] <- 'Gochnatioideae'
    }
  else if (df_sum$Tribu[i] %in% mutisioideae) {
    df_sum$subfam[i] <- 'Mutisioideae'
    }
  else if (df_sum$Tribu[i] %in% stifftioideae) {
    df_sum$subfam[i] <- 'Stifftioideae'
    }
  else{
    df_sum$subfam[i] <- 'Barnadesioideae'
  }
}

df_sum$subfam <- factor(as.factor(df_sum$subfam), levels = c( 'Cichorioideae', 'Mutisioideae','Gochnatioideae','Stifftioideae', 'Asteroideae', 'Barnadesioideae'))
p5 <- ggplot(data = df_sum, aes (y = sv, x = subfam, color = subfam))
p5 <- p5 + geom_boxplot() + scale_color_viridis_d(option = 'A', end = 0.9) + scale_x_discrete(labels = NULL) + ggtitle('SV ~ Subfamilia') + theme_bw()
class(df_sum$Latitud)
p4 + p5
Latitud <- gsub(',','.', df_sum$Latitud)
df_sum$Latitud <- as.numeric(Latitud)
Longitud <- gsub(',','.', df_sum$Longitud)
df_sum$Longitud <-  as.numeric(Longitud)
head(df_sum)



p6 <- ggplot(data = df_sum, aes(y = sv, x = Latitud))
p6 + geom_point() + geom_smooth()
df_sum_subdf <- filter(df_sum, muestreo == 'norte_2022' | muestreo == 'norte_2020')
p6 <- ggplot(data = df_sum_subdf, aes(y = sv, x = Latitud))
p6 + geom_point() + geom_smooth()
df_sum_subdf <- filter(df_sum, muestreo == 'norte_2022')
p6 <- ggplot(data = df_sum_subdf, aes(y = sv, x = Latitud))
p6 + geom_point() + geom_smooth()
unique(df_sum_subdf$cod_gps)
df_sum_subdf <- filter(df_sum_subdf, cod_gps != 'Pisco Huasi' )

p7 <- ggplot(data = df_sum_subdf, aes(y = sv, x = Latitud))
p7 + geom_point() + geom_smooth()

### mapas
### https://kdmurray.id.au/post/2021-05-06_how-to-extract-points-from-raster-r/

library(raster)

### en realidad son el mismo punto
df_sum[58,'Latitud'] = -50.4996296
df_sum[58,'Longitud'] = -73.0435309
df_sum[59,'Latitud'] = -50.4925138
df_sum[59,'Longitud'] = -73.0476839

#setwd('./tablas segundo capitulo')
getwd()
load_bioclim = function(i) {
  raster(sprintf("wc2.1_30s_bio/wc2.1_30s_bio_%02d.tif", i),
         varname=sprintf("Bio%02d", i))
}

rasters = do.call("stack", lapply(1:19, load_bioclim))
coord = df_sum
coordinates(coord) = ~ Longitud + Latitud

point.values = raster::extract(rasters, coord)
class(df_sum)
df_sum$bio1 <- point.values[,1]
df_sum$bio4 <- point.values[,4]
df_sum$bio12 <- point.values[,12]
df_sum$bio15 <- point.values[,15]






p_bio1 <- ggplot(data = df_sum, aes(y = sv, x = bio1, color = muestreo))
p_bio1 + geom_point() + geom_smooth() + facet_wrap(~ muestreo) + ggtitle('SV ~ temperatura') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

p_bio1 <- ggplot(data = df_sum, aes(y = sv, x = bio1))
p_bio1 + geom_point() + geom_smooth() + ggtitle('SV ~ temperatura') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

p_bio4 <- ggplot(data = df_sum, aes(y = sv, x = bio4, color = muestreo))
p_bio4 + geom_point() + geom_smooth() + facet_wrap(~ muestreo) + ggtitle('SV ~ variabilidad en temp') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()


p_bio4 <- ggplot(data = df_sum, aes(y = sv, x = bio4))
p_bio4 + geom_point() + geom_smooth() + ggtitle('SV ~ variabilidad en temp') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

       
p_bio12 <- ggplot(data = df_sum[df_sum$bio12 < 1300,], aes(y = sv, x = bio12))
p_bio12 + geom_point() + geom_smooth() + ggtitle('SV ~ precipitaciones') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()


p_bio15 <- ggplot(data = df_sum, aes(y = sv, x = bio15))
p_bio15 + geom_point() + geom_smooth()  + ggtitle('SV ~ var precipitaciones') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

## cobertura de arboles
# unzip files 
#for (i in 1:2) {
#  for (z in 1:6) {
#    unzip(sprintf("cover/gm_ve_v2_%1d_%1d.zip", i,z), exdir = 'cover')
#  }
#}
## Esto fue para armar la imagen cover.tif ahora solo la cargo
#lista = list()
#for (i in 1:2) {
#  for (z in 1:6) {
#    if (i == 1) {
#    lista [[z]] = raster(sprintf("cover/gm_ve_v2_%1d_%1d.tif", i, z),
#                         varname=sprintf("cover1-%1d", z))
#    }
#    else {
#      lista [[z+6]] = raster(sprintf("cover/gm_ve_v2_%1d_%1d.tif", i,z),
#                           varname=sprintf("cover2-%1d", z))
#    }
#  }
#}
#lista$filename = 'cover.tif'
#lista$overwrite = TRUE
#m <- do.call(raster::merge, lista)
getwd()
m = raster('cover.tif', varname = 'tree.cover')

coord = df_sum
coordinates(coord) = ~ Longitud + Latitud
point.values = raster::extract(m, coord)

df_sum$cover = point.values

p_cover <- ggplot(data = df_sum, aes(y = sv, x = cover))
p_cover + geom_point() + geom_smooth() + ggtitle('SV ~ tree cover') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

df_sum2 = df_sum %>% filter(cover <= 15)
p_cover2 <- ggplot(data = df_sum2, aes(y = sv, x = cover))
p_cover2 + geom_point() + geom_smooth() + ggtitle('SV ~ tree cover') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

### wind speed
#10
ws = raster('ARG_wind-speed_10m.tif', varname = 'wind_speed')
point.values = raster::extract(ws, coord)
table(as.factor(point.values))
df_sum$wind.speed = point.values
p_ws <- ggplot(data = df_sum, aes(y = sv, x = wind.speed))
p_ws + geom_point() + geom_smooth() + ggtitle('SV ~ wind speed') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

# 50
ws_50 = raster('ARG_wind-speed_50m.tif', varname = 'wind_speed')
point.values = raster::extract(ws_50, coord)
df_sum$wind.speed50 = point.values
p_ws_50 <- ggplot(data = df_sum, aes(y = sv, x = wind.speed50))
p_ws_50 + geom_point() + geom_smooth() + ggtitle('SV ~ wind speed50m') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

p_ws_50 <- ggplot(data = df_sum, aes(y = l_papus, x = wind.speed50))
p_ws_50 + geom_point() + geom_smooth() + ggtitle('l_papus ~ wind speed50m') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()



### air density
ad = raster('ARG_air-density_10m.tif', varname = 'air_density')
point.values = raster::extract(ad, coord)
df_sum$air.density = point.values
p_ad <- ggplot(data = df_sum, aes(y = l_papus, x = air.density))
p_ad + geom_point() + geom_smooth() + ggtitle('SV ~ air.density') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
p_ad <- ggplot(data = df_sum, aes(y = sv, x = air.density))
p_ad + geom_point() + geom_smooth() + ggtitle('SV ~ air.density') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()


### capacity factor IEC 
iec = raster('ARG_cf_IEC1.tif', varname = 'capacity_factor_IEC')
point.values = raster::extract(iec, coord)
df_sum$cf_iec = point.values
p_iec <- ggplot(data = df_sum, aes(y = sv, x = cf_iec))
p_iec + geom_point() + geom_smooth() + ggtitle('SV ~ iec') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()


### parto el set de datos por las tribus que mas representadas estan

stribus = names(which(table(df_sum$Tribu) > 6))
subdf_tribes = df_sum[df_sum$Tribu %in% stribus,]

#por ecoregion
count
c_ecoregion <- subdf_tribes %>% group_by(Ecoregion) %>% summarise(mSV=max(sv), n = n())
p_tribes <- ggplot(data = subdf_tribes , aes (y = sv, x = Ecoregion, color = Ecoregion))
p_tribes + geom_boxplot() + geom_text(data = c_ecoregion, aes(x= Ecoregion, y = 5e-08+mSV, label = n ), show.legend = F) + scale_color_viridis_d(option = 'A', end = 0.90) + ggtitle('SV ~ Ecoregion') + theme_bw()

 
p_tribes_cover <- ggplot(data = subdf_tribes, aes(y = sv, x = cover))
p_tribes_cover + geom_point() + geom_smooth() + ggtitle('SV ~ tree cover') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
# wind speed
p_ws <- ggplot(data = subdf_tribes, aes(y = sv, x = wind.speed))
p_ws + geom_point() + geom_smooth() + ggtitle('SV ~ wind speed') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
# air density
p_ad <- ggplot(data =subdf_tribes, aes(y = sv, x = air.density))
p_ad + geom_point() + geom_smooth() + ggtitle('SV ~ air.density') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
### Por subfamilia

subfamilies = names(which(table(df_sum$subfam) > 10))
subdf_sf = df_sum[df_sum$subfam %in% subfamilies,]

c_ecoregion <- subdf_sf %>% group_by(Ecoregion) %>% summarise(mSV=max(sv), n = n())
p_sf <- ggplot(data = subdf_sf , aes (y = sv, x = Ecoregion, color = Ecoregion))
p_sf + geom_boxplot() + geom_text(data = c_ecoregion, aes(x= Ecoregion, y = 2e-08+mSV, label = n )) + scale_color_viridis_d(option = 'A', end = 0.90) + ggtitle('SV ~ Ecoregion') + theme_bw()

## cover
p_sf_cover <- ggplot(data = subdf_sf, aes(y = sv, x = cover))
p_sf_cover + geom_point() + geom_smooth() + ggtitle('SV ~ tree cover') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

## wind speed
p_ws <- ggplot(data = subdf_sf, aes(y = sv, x = wind.speed))
p_ws + geom_point() + geom_smooth(method = 'lm') + ggtitle('SV ~ wind speed') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
## air density

p_ad <- ggplot(data =subdf_sf, aes(y = sv, x = air.density))
p_ad + geom_point() + geom_smooth() + ggtitle('SV ~ air.density') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

### export datframe final
#write.csv(df_sum, 'final_df.csv')

## solo con setosas ################# 
df_sum2 = df_sum %>% filter(tipo_papus == 'setose')


### temperatura #############
p_bio1 <- ggplot(data = df_sum2, aes(y = sv, x = bio1, color = muestreo))
p_bio1 + geom_point() + geom_smooth() + facet_wrap(~ muestreo) + ggtitle('SV ~ temperatura') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

p_bio1 <- ggplot(data = df_sum2, aes(y = sv, x = bio1))
p_bio1 + geom_point() + geom_smooth() + ggtitle('SV ~ temperatura') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
### var temperatura ####################
p_bio4 <- ggplot(data = df_sum2, aes(y = sv, x = bio4, color = muestreo))
p_bio4 + geom_point() + geom_smooth() + facet_wrap(~ muestreo) + ggtitle('SV ~ variabilidad en temp') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

p_bio4 <- ggplot(data = df_sum2, aes(y = sv, x = bio4))
p_bio4 + geom_point() + geom_smooth() + ggtitle('SV ~ variabilidad en temp') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

### precipitaciones ########
p_bio12 <- ggplot(data = df_sum2[df_sum2$bio12 < 1300,], aes(y = sv, x = bio12))
p_bio12 + geom_point() + geom_smooth() + ggtitle('SV ~ precipitaciones') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

### var precipitaciones ###################
p_bio15 <- ggplot(data = df_sum2, aes(y = sv, x = bio15))
p_bio15 + geom_point() + geom_smooth()  + ggtitle('SV ~ var precipitaciones') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

### cover #########
m = raster('cover.tif', varname = 'tree.cover')
coord = df_sum2
coordinates(coord) = ~ Longitud + Latitud
point.values = raster::extract(m, coord)

df_sum2$cover = point.values

p_cover <- ggplot(data = df_sum2, aes(y = sv, x = cover))
p_cover + geom_point() + geom_smooth() + ggtitle('SV ~ tree cover') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

### wind speed #######
#10
ws = raster('ARG_wind-speed_10m.tif', varname = 'wind_speed')
point.values = raster::extract(ws, coord)
table(as.factor(point.values))
df_sum2$wind.speed = point.values
p_ws <- ggplot(data = df_sum2, aes(y = sv, x = wind.speed))
p_ws + geom_point() + geom_smooth() + ggtitle('SV ~ wind speed') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

# 50
ws_50 = raster('ARG_wind-speed_50m.tif', varname = 'wind_speed')
point.values = raster::extract(ws_50, coord)
df_sum2$wind.speed50 = point.values
p_ws_50 <- ggplot(data = df_sum2, aes(y = sv, x = wind.speed50))
p_ws_50 + geom_point() + geom_smooth() + ggtitle('SV ~ wind speed50m') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

p_ws_50 <- ggplot(data = df_sum2, aes(y = l_papus, x = wind.speed50))
p_ws_50 + geom_point() + geom_smooth() + ggtitle('l_papus ~ wind speed50m') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()



### air density ####### 
ad = raster('ARG_air-density_10m.tif', varname = 'air_density')
point.values = raster::extract(ad, coord)
df_sum2$air.density = point.values
p_ad <- ggplot(data = df_sum2, aes(y = l_papus, x = air.density))
p_ad + geom_point() + geom_smooth() + ggtitle('SV ~ air.density') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()
p_ad <- ggplot(data = df_sum2, aes(y = sv, x = air.density))
p_ad + geom_point() + geom_smooth() + ggtitle('SV ~ air.density') +
  scale_color_viridis_d(option = 'A', direction = -1)  + theme_bw()

## mapas ################## 
setwd('./resumen')

library(sf);library(spData);library(spDataLarge)
library(tmap);library(leaflet); library(patchwork)

argentina <- st_read('provincias.shp')
argentina <- filter(argentina, FNA != "Provincia de Tierra del Fuego, Antártida e Islas del Atlántico Sur")
head(df_sum)
points <- st_as_sf(df_sum, 
                   coords = c( x = 'Longitud', y = 'Latitud'),
                   crs = 'WGS84'
)
mapa.arg <- tm_shape(argentina) + tm_polygons()
mapa.points <- mapa.arg + tm_shape(points) + tm_dots(size = 1, col = 'Ecoregion')
mapa.points

df_sum$cod_gps = as.factor(df_sum$cod_gps)
df_sum$Ecoregion = as.factor(df_sum$Ecoregion)
df_sum$muestreo = as.factor(df_sum$muestreo)
subdf_sum2 <- data.frame(sitio = factor(rep('', length(levels(df_sum$cod_gps))), levels = levels(df_sum$cod_gps)) ,
                     ecoregion = factor(rep('', length(levels(df_sum$cod_gps))), levels = levels(df_sum$Ecoregion)),
                     long = rep(0, length(levels(df_sum$cod_gps))),
                     lat = rep(0, length(levels(df_sum$cod_gps))),
                     alt = rep(0, length(levels(df_sum$cod_gps))),
                     n = rep(0, length(levels(df_sum$cod_gps))),
                     muestreo = factor(rep('', length(levels(df_sum$cod_gps))), levels = levels(df_sum$muestreo))
)

#df_sum[df_sum$Tribu.Asteraceae == 'Astereae',][1,4] # prueba
for (i in 1:length(levels(df_sum$cod_gps))){
  subdf_sum2$sitio[i] <- levels(df_sum$cod_gps)[i]
  subdf_sum2$ecoregion[i] <- df_sum[df_sum$cod_gps == subdf_sum2$sitio[i],]$Ecoregion[1]
  subdf_sum2$long[i] <- df_sum[df_sum$cod_gps == subdf_sum2$sitio[i],]$Longitud[1] # como todo el sitio tiene misma long, lat y alt selecciono el priemero que salga
  subdf_sum2$lat[i] <- df_sum[df_sum$cod_gps == subdf_sum2$sitio[i],]$Latitud[1]
  subdf_sum2$alt[i] <- df_sum[df_sum$cod_gps == subdf_sum2$sitio[i],]$altitud[1]
  subdf_sum2$n[i] <- length(which(df_sum$cod_gps == subdf_sum2$sitio[i]))
  subdf_sum2$muestreo[i] <- df_sum[df_sum$cod_gps == subdf_sum2$sitio[i],]$muestreo[1]
}
df_lau = subdf_sum2[,c('sitio','long','lat','ecoregion', 'n')]
write.csv(df_lau,'mapa.csv')
points_sitios <- st_as_sf(subdf_sum2, 
                          coords = c( x = 'long', y = 'lat'),
                          crs = 'WGS84'
)
mapa.arg <- tm_shape(argentina) + tm_polygons()
mapa.points <- mapa.arg + tm_shape(points_sitios) + 
  tm_symbols(col = 'ecoregion',size = 'n', alpha = 0.5, 
             palette = viridis(3, directio = -1, begin = 0.2), 
             title.col = 'Ecoregion', title.size = 'Riqueza') 
 
#pdf('mapa-muestreo.pdf')
setwd('..')

png('mapa-muestreo_act.png', width = 14, height = 17, units = 'cm', res = 300)
mapa.points  +  tm_graticules(alpha = 0.3, lwd = 0.5) + tm_compass(position = c("left", "top"), size = 2) + 
  tm_scale_bar(breaks = c(0, 200, 400, 600), text.size = 0.5, position = c('right', 'bottom')) +
  tm_layout(outer.margins = c(0.001,0.001,0.001,0.001),
            legend.title.size = 0.7,
            legend.text.size = 0.5,
            legend.outside = T)
dev.off()


### independent effects ##################
library(hier.part)
colnames(df_sum)

corr <-cor(na.omit(df_sum[,c('masa', 'area_p', 'pl_m','cover', 'wind.speed', 'wind.speed50', 'air.density', 'bio1', 'bio12','Latitud', 'Longitud')], 1))

corr
library(viridis)
ggcorrplot(corr, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T)
## tribus muestreadas 
df_sum = df_sum[!df_sum$Tribu == 'Coreopsideae',] ## saco thelesperma
l <- reorder(levels(df_sum$Tribu), as.vector(table(df_sum$Tribu)), decreasing = T) # creo un factor con los levels ordenados
levels(l)
df_sum$Tribu <- factor(df_sum$Tribu, levels = levels(l))
my_theme <-   theme(
  legend.text = element_text(size = 8),
  legend.background = element_rect(fill = "white"),
  plot.title = element_text(vjust = -1, hjust = 0.5),
  plot.subtitle = element_text(vjust = -1, hjust = 0.5),
  axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
  axis.text.y = element_text(size = 8),
  axis.title = element_text(size = 10),
  axis.line = element_line(colour = "grey50"),
  plot.margin = margin(6,0,0,6),
  panel.grid = element_line(color = "#b4aea9"),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(linetype = "dashed"),
  panel.background = element_rect(fill = 'white', color = "white"),
  plot.background = element_rect(fill = "white", color = "white")
)
bar.plot <- ggplot(df_sum, aes(x = Tribu))
tribes.dist <- bar.plot + geom_bar(aes(fill = Tribu), show.legend = F) + 
  scale_y_continuous(expand = c(0,0), breaks = seq(0,15, by = 3), limits = c(0,15)) + 
  scale_fill_viridis_d('Tribus', option = 'A', direction = -1) + xlab('Tribus') + ylab('n especies') + my_theme
#png('tribus.png', width = 15, height = 10, units = 'cm', res = 300)
tribes.dist
dev.off()


eco_count <- count(df_sum, Ecoregion)
df_sum$Ecoregion <- factor(df_sum$Ecoregion, 
                       levels = eco_count$Ecoregion[order(eco_count$n, decreasing = T)])

bar.plot <- ggplot(df_sum, aes(x = Ecoregion, fill = Tribu))

#png('ecoreg_trib.png', width = 18, height = 18, units = 'cm', res = 300)
bar.plot + geom_bar() + scale_y_continuous(expand = c(0,0), 
                                           breaks = seq(0,50, by = 10),
                                           limits = c(0,50)) + 
  scale_fill_viridis_d('Tribus', option = 'A', direction = -1) + xlab('Ecoregiones') +ylab('n especies') + 
  my_theme + theme(axis.text.x = element_text(size = 8, angle = 15, vjust = 0.6))
dev.off()






df_sum_ord = df_sum[,c('wind.speed', 'air.density', 'bio1', 'bio12','Latitud', 'Longitud', 'Ecoregion', 'Tribu','tipo_papus','masa', 'area_p')]
hier.part(df_sum$sv, df_sum_ord, fam = "gaussian", link = "identity", gof = 'logLik', plot = F)

df_py = df_sum[,c('sv','wind.speed', 'bio1', 'bio12','Latitud', 'Longitud', 'Ecoregion', 'Tribu', 'habito','tipo_papus','masa', 'area_p', 'pl_m')]
write_csv(df_py,'df_ev.csv' )


### graficos de sv en funcion de ciertas variables


max(df_sum$sv)
min(df_sum$sv)
sd(df_sum$sv)
mean(df_sum$sv)
p_sv2 <- ggplot(df_sum, aes(x = sv))
p_sv2 <- p_sv2 + geom_histogram(aes(x = sv, y = after_stat(density)),
                        col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                        fill = 'white', bins = 30) + scale_x_continuous(breaks = seq(0,3.5, by = 0.5)) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) +
  xlab('Velocidad de Asentamiento (m/s)') + ylab('Densidad') +
  scale_y_continuous(expand = c(0,0))+ 
  theme(
    axis.text.x = element_text(size = 12, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(colour = "grey50"),
    plot.margin = margin(6,0,0,6),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = 'white', color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  )


p_sv2

dist_sv = wrap_plots(p_sv + ggtitle('a)'), p_sv2 + ggtitle('b)')) 
pdf('dist_sv.pdf')
dist_sv
dev.off()
png('dist_sv.png',width = 22, height = 13, units = 'cm', res = 600)
dist_sv
dev.off()

p = ggplot(df_sum, aes(x = area_p, y = sv, color = tipo_papus))
p + geom_point() + geom_smooth(method = 'lm')

p = ggplot(df_sum, aes(x = log(masa), y = sv, color = tipo_papus))
p + geom_point() + geom_smooth()

df_sum$pl_m2 = df_sum$masa/df_sum$area_p

p = ggplot(df_sum, aes(x = sqrt(pl_m2), y = sv, color = tipo_papus))
p + geom_point() + geom_smooth()

p = ggplot(df_sum, aes(x = sqrt(pl_m), y = masa, color = tipo_papus))
p + geom_point() + geom_smooth()

p = ggplot(df_sum, aes(x = sqrt(pl_m), y = area_p, color = tipo_papus))
p + geom_point() + geom_smooth()

p = ggplot(df_sum, aes(x = masa, y = area_p, color = tipo_papus))
p + geom_point() + geom_smooth(method = 'lm')


### modelos #######


library(lme4)

model1 = lmer(sv ~ tipo_papus + (1|Tribu), data = df_sum)
summary(model1)
model1 = lmer(sv ~ tipo_papus + (1|Tribu), data = df_sum)
model0 = lmer(sv ~ (1|Tribu), data = df_sum)
anova(model1,model0)
install.packages("emmeans")
#install.packages("multcomp")
library(emmeans)
library(multcomp)
emmeans_model <- emmeans(model1, specs = "tipo_papus")
posthoc_results <- pairs(emmeans_model, adjust = "tukey")
summary(posthoc_results)

mean(df_sum[df_sum$tipo_papus == 'setose',]$sv)
mean(df_sum[df_sum$tipo_papus == 'paleaceous',]$sv)
mean(df_sum[df_sum$tipo_papus == 'epappose',]$sv)
mean(df_sum[df_sum$tipo_papus == 'aristate',]$sv)




letras <- c("a", "a", "b", "c")
names(letras) <- c("aristate", "epappose", "paleaceous", "setose")
## los mismos colores que el paper
colores <- c("aristate" = "#AA337DFF", 
             "epappose" = "#000004FF", 
             "paleaceous" = "#F7725CFF", 
             "setose" = "#FDE2A2FF")

n_counts <- df_sum %>%
  group_by(tipo_papus) %>%
  summarise(n = n())

p_boxplot = ggplot(df_sum, aes(x = tipo_papus, y = sv, fill = tipo_papus)) +
  geom_boxplot() +
  scale_fill_manual(values = colores) + 
  geom_text(data = data.frame(tipo_papus = names(letras), label = letras), 
            aes(x = tipo_papus, y =c(3.03, 3.03, 2.5, 2.3), label = label), 
            size = 3, vjust = 0) +  
  geom_text(data = n_counts, 
            aes(x = tipo_papus, y = min(df_sum$sv) - 1e-02, label = paste0("n=", n)), 
             size = 3, vjust = 1.5, color = "black") + 
  labs(x = "Tipo de Papus", y = "VA (m/s)") +
  my_theme + theme(axis.text.x = element_text(angle = 0),legend.position = "none") + 
  ggtitle("(a)") + theme(plot.title = element_text(hjust = -0.05,vjust = 0.5))

#png('tipo_papus_boxplot.png', width = 15, height = 10, units = 'cm', res = 300)
p_boxplot
dev.off()

model1 = lmer(sv ~ tipo_papus_b + (1|Tribu), data = df_sum)
summary(model1)
model0 = lmer(sv ~ (1|Tribu), data = df_sum)
anova(model1,model0)

mean(df_sum[df_sum$tipo_papus_b == 'setose',]$sv)
mean(df_sum[df_sum$tipo_papus_b == 'other',]$sv)


model1 = lmer(sv~ scale(area_p) + scale(masa) + tipo_papus - 1 + (1|Tribu), data = df_sum)
summary(model1)
model1b = lmer(log(sv)~ scale(masa) + tipo_papus + (1|Tribu), data = df_sum)
model1c = lmer(log(sv)~ scale(masa) + scale(area_p) + (1|Tribu), data = df_sum)

summary(model1)
anova(model1,model1b) ## area no significativa
anova(model1,model1c) ## tipo papus significtivo

### SV en funcion de variables ambientales y tipo de papus 

df_sum$tipo_papus_b = ifelse(df_sum$tipo_papus == 'setose', 'setose', 'other')
summary(as.factor(df_sum$tipo_papus_b))


papus = df_sum %>% select(sv,tipo_papus) %>% group_by(tipo_papus) %>% summarise_all(list(mean,sd)) 
papusb = df_sum %>% select(sv,tipo_papus_b) %>% group_by(tipo_papus_b) %>% summarise_all(list(mean,sd)) 



p_boxplot2 = ggplot(df_sum, aes(x = tipo_papus_b, y = sv, fill= tipo_papus_b)) +
  geom_boxplot() + scale_fill_manual(values = magma(2,direction = 1, end = 0.93)) + 
  labs(x = "Tipo de Papus", y = "VA (m/s)") +
  my_theme + theme(axis.text.x = element_text(angle = 0),legend.position = "none") + ggtitle("(b)") +
  theme(plot.title = element_text(hjust = -0.05,vjust = 0.5))


design <- "AAABB"
wrap_plots(A = p_boxplot, B = p_boxplot2, design = design) 




p_ws <- ggplot(data = df_sum, aes(y = sv, x = wind.speed, color = tipo_papus_b))
p_ws + geom_point() + geom_smooth(method = 'lm')  +
  scale_color_viridis_d(option = 'A', direction = 1, end = 0.95)  + theme_bw()

library(lme4)
model1 = lmer(sv ~  wind.speed*tipo_papus_b + (1|Tribu), data = df_sum)
summary(model1)
model1.1 = lmer(sv ~  wind.speed + tipo_papus_b + (1|Tribu), data = df_sum)
anova(model1, model1.1)
model1.2 = lmer(sv ~  tipo_papus_b + (1|Tribu), data = df_sum)
anova(model1.1, model1.2) ## es lo mismo simplificando el modelo. Como sabiamos tambien por el random forest.

df_sum2 = df_sum[df_sum$tipo_papus_b == 'setose',] 

corr <-cor(df_sum2[,c('masa', 'area_p', 'pl_m')])


df_sum2[which(df_sum2$masa == max(df_sum2$masa)),]
df_sum2[which(df_sum2$area_p == max(df_sum2$area_p)),]

df_sum2 = df_sum2[!df_sum2$spp == 'Thelesperma_megapotamicum',] 

### independent effects 
library(hier.part)
hier.part(df_sum2$sv, df_sum2[,c('bio1', 'bio4', 'cover', 'wind.speed')], fam = "gaussian", link = "identity", gof = 'logLik', plot = T)

plot(log(sv)~sqrt(pl_m), data = df_sum2)
df_sum2$sqrt_pl_m = sqrt(df_sum2$pl_m)
library(performance)
model01 = lmer(sv ~ scale(area_p)  + scale(masa) + (1|Tribu), data = df_sum2)
summary(model01)
modelo01.1 = lmer(sv ~ scale(area_p)   + (1|Tribu), data = df_sum2)
anova(model01, modelo01.1)
modelo01.2 = lmer(sv ~  scale(masa) + (1|Tribu), data = df_sum2)
anova(model01, modelo01.2)


library(ggeffects)
preds = ggpredict(model01, c('area_p[all]'))
p_ap = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "area del papus (mm2)", y = "") + geom_point(data = df_sum2, aes(x = area_p, y = sv)) + theme_bw()
p_ap = p_ap + annotate('text',x = 700, y = 2.2, label = 'p > 0.05')

preds = ggpredict(model01, c('masa[all]'))
p_m = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "masa (mg)", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = masa, y = sv)) + theme_bw()
p_m = p_m + annotate('text',x = 7.5, y = 2.2, label = 'p < 0.05 ***')


r2_values <- r2(model01)
print(r2_values)

model1 = lmer(sv ~ scale(sqrt_pl_m) + (1|Tribu), data = df_sum2)
summary(model1)
model0 = lmer(sv ~(1|Tribu), data = df_sum2)
anova(model1, model0)
r2(model1)
preds = ggpredict(model1, c('sqrt_pl_m[all]'))
p_pl = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "raiz cuadrada de carga de la pluma", y = "") + geom_point(data = df_sum2, aes(x = sqrt_pl_m, y = sv)) + theme_bw()
p_pl = p_pl + annotate('text',x = 0.15, y = 2.4, label = 'p < 0.05 ***')

wrap_plots(p_m + ggtitle('(a)'),p_ap +  ggtitle('(b)'), p_pl + ggtitle('(c)'))
library(arm)
retribu <- ranef(model1)
setribu <- se.ranef(model1)
raneftribu = data.frame(
  int = retribu$Tribu[,1],
  se = setribu$Tribu[,1]
)
## para graficar los efectos aleatorios
p_r = ggplot(raneftribu,aes(x= rownames(raneftribu),y=int)) +
  geom_point(size=2.5, col="red") +
  geom_errorbar(aes(ymin=int-se, ymax=int+se),width=.5)+ geom_hline(yintercept = 0, color = "blue") +
  theme_bw() +labs(y="Intercepto", x="Tribu")+coord_flip()
library(patchwork)


### modelos de caracteres extr[insecos]
colnames(df_sum2)
model1 = lmer(sv ~  scale(wind.speed) + scale(bio12) + scale(cover) + (1|Tribu), data = df_sum2)
summary(model1)
model1.1 = lmer(sv ~ scale(wind.speed) + scale(cover)  + (1|Tribu), data = df_sum2)
anova(model1, model1.1)
model1.2 = lmer(sv ~ scale(bio12) + scale(cover) + (1|Tribu), data = df_sum2)
anova(model1, model1.2)
model1.3 = lmer(sv ~ scale(bio12) + scale(wind.speed) + (1|Tribu), data = df_sum2)
anova(model1, model1.3)
r2_values <- r2(model1)
r2_values
library(ggeffects)

preds = ggpredict(model1, 'wind.speed[all]')
p_w = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Velocidad del viento (m/s)", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = wind.speed, y = sv)) + theme_bw()
p_w = p_w + annotate('text',x = 3, y = 2.2, label = 'p < 0.05 *') + ggtitle('(b)')


preds = ggpredict(model1, 'bio12[all]')
p_p = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Precipitación anual media (mm)", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = bio12, y = sv)) + theme_bw()
p_p = p_p + annotate('text',x = 750, y = 2.2, label = 'p < 0.05 *') + ggtitle('(a)')

preds = ggpredict(model1, 'cover[all]')
p_c = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Cobertura de árboles", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = cover, y = sv)) + theme_bw()
p_c = p_c + annotate('text',x = 40, y = 2.2, label = 'p > 0.05') + ggtitle('(c)')

library(arm)
retribu <- ranef(model1)
setribu <- se.ranef(model1)
raneftribu = data.frame(
  int = retribu$Tribu[,1],
  se = setribu$Tribu[,1]
)
## para graficar los efectos aleatorios
p_r = ggplot(raneftribu,aes(x= rownames(raneftribu),y=int)) +
  geom_point(size=2.5, col="red") +
  geom_errorbar(aes(ymin=int-se, ymax=int+se),width=.5)+ geom_hline(yintercept = 0, color = "blue") +
  theme_bw() +labs(y="Intercepto", x="Tribu")+coord_flip() + ggtitle('(d)')

(p_p + p_w)/(p_c + p_r)




library(nlme)
modelo <- lme(fixed = log(sv) ~ wind.speed, data = df_sum2,
              random = ~ 1 | Tribu,
              correlation = corCompSymm(form = ~ 1 | Tribu))
summary(modelo)
plot(modelo)


model2 = lmer(log(sv) ~   wind.speed + (1|Tribu), data = df_sum2)
preds = ggpredict(model2, 'wind.speed[all]')
p_f = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Velocidad del viento (m/s)", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = wind.speed, y = sv)) + theme_bw()
library(arm)
retribu <- ranef(model1)
setribu <- se.ranef(model1)
raneftribu = data.frame(
  int = retribu$Tribu[,1],
  se = setribu$Tribu[,1]
)

p_r = ggplot(raneftribu,aes(x= rownames(raneftribu),y=int)) +
  geom_point(size=2.5, col="red") +
  geom_errorbar(aes(ymin=int-se, ymax=int+se),width=.5)+ geom_hline(yintercept = 0, color = "blue") +
  theme_bw() +labs(y="Intercepto", x="Tribu")+coord_flip()
library(patchwork)
p_f + p_r


df_sum2 = df_sum2[!df_sum2$masa == max(df_sum2$masa),]
df_sum2 = df_sum2[!df_sum2$area_p == max(df_sum2$area_p),]
model2 = lmer(log(sv) ~   masa + l_papus + bio12 + wind.speed + (1|Tribu), data = df_sum2)
summary(model2)

preds = ggpredict(model2, 'masa[all]')
p_f = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "masa", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = wind.speed, y = sv)) + theme_bw()
library(arm)
retribu <- ranef(model2)
setribu <- se.ranef(model2)
raneftribu = data.frame(
  int = retribu$Tribu[,1],
  se = setribu$Tribu[,1]
)

p_r = ggplot(raneftribu,aes(x= rownames(raneftribu),y=int)) +
  geom_point(size=2.5, col="red") +
  geom_errorbar(aes(ymin=int-se, ymax=int+se),width=.5)+ geom_hline(yintercept = 0, color = "blue") +
  theme_bw() +labs(y="Intercepto", x="Tribu")+coord_flip()
library(patchwork)
p_f + p_r

preds = ggpredict(model2, 'l_papus[all]')
p_f = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "largo del papus", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = wind.speed, y = sv)) + theme_bw()

p_f
