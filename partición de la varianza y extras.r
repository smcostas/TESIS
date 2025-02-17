library(tidyverse);library(viridis);library(lme4); library(fitdistrplus); library(MuMIn); library(cAIC4); library(parallel);library(pbkrtest)
library(ggcorrplot); library(PerformanceAnalytics); library(reshape2); library(hier.part); library(patchwork)
 ## Tabla ####
df <- read.csv('datos_finales.csv', header = T)
df$tiempo <- df$tiempo/1000000
df$sv <- 2.2/df$tiempo
plot(log(df$sv)~log(df$masa))

df2 <- df

df2$area_p<-df2$l_papus^2*pi
df2$vol_c<-pi*(df2$a_cisela/2)^2*(df2$l_cipsela/2)

df2$pl_m <- df2$masa/df2$area_p
df2$pl_v <- df2$vol_c/df2$area_p

head(df2)

### distribución de la masa ####
p <- ggplot(df2, aes(x = masa)) + geom_histogram(aes(x = masa),
                        col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                        fill = 'white', bins = 35) +scale_x_continuous(breaks=seq(0,3, by = 0.25))  + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) + scale_x_continuous(breaks=seq(0,3, by = 0.25)) + 
  xlab('masa') + ylab('Densidad')
  scale_y_continuous(expand = c(0,0))+ my_theme
p









head(df2)
df_sum <- df2 %>% dplyr::select(pob, indiv,fruto, sv, pl_m, pl_v, masa, l_cipsela, a_cisela, l_papus, area_p, vol_c, capitulos, altura_planta ) %>% ## selecciono varaibles de interes
  group_by(pob,indiv)  %>% summarise_all(lst(mean, sd))# marcamos las variables que agrupan
df_sum
corr <- round(cor(na.omit(df_sum[,c(10,11,8,9,7,12,5,6,13,14)])), 1) ## tabla de correlaciones
df_sum[,12]
my_theme <- theme(
  axis.text.x = element_text(size = 12, vjust = 0.5),
  axis.text.y = element_text(size = 12, vjust = 0.5),
  axis.line = element_line(colour = "grey50"),
  plot.margin = margin(6,0,0,6),
  panel.grid = element_line(color = "#b4aea9"),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(linetype = "dashed"),
  panel.background = element_rect(fill = 'white', color = "white"),
  plot.background = element_rect(fill = "white", color = "white"),
  legend.title=element_blank()
)

## Correlacion ######
#no hay autocorrelacion entre las variables (sacando las variables que fueron creadas)
ggcorrplot(corr, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T)



## distribución del Settling velocity ####
head(df_sum)
p <- ggplot(df_sum, aes(x = sv_mean))
p <- p + geom_histogram(aes(x = sv_mean, y = after_stat(density)),
                        col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                        fill = 'white', bins = 30) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) +
  xlab('Velocidad de Asentamiento media por fruto (m/s)') + ylab('Densidad') +
  scale_y_continuous(expand = c(0,0))+ my_theme
p
p2 <- ggplot(df_sum, aes(x = sv_sd)) ### error
p2 <- p2 + geom_histogram(aes(x = sv_sd, y = ..density..),
   col = viridis(1, begin = 0.3, option = 'A', direction = 1),
  fill = 'white', bins = 30) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) +
  xlab('sd de Velocidad de Asentamiento (m/s)') + ylab('Densidad') +
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
p2

mean(print(df_sum$sv_sd/df_sum$sv_mean)) ### el coeficiente de variación por fruto es bajisimo. habla de que en cada tirada es practicamente igual 

### me pregunto si no habría que hacer lo mismo para entre frutos dentro de individuos y para individuos dentro de poblaciones ahora lo vamos a ver en particion de la varianza

## modelos mixtos para ver la particion de la varianza ####
df_sum_f <- df2 %>% dplyr::select(pob, indiv, fruto, sv, pl_m, pl_v, masa, l_cipsela, a_cisela, l_papus, area_p, vol_c, capitulos, altura_planta ) %>% ## selecciono varaibles de interes
  group_by(pob,indiv,fruto)  %>% summarise_all(lst(mean))# marcamos las variables que agrupan
df_sum_f
### SV ####
m1_sv = lmer(sv_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = T)
summary(m1_sv) ## la residual incluyo intra y entre frutos
### PL_m ####
m2_pl_m = lmer(pl_m_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)
summary(m2_pl_m)
### PL_v ####
m3_pl_v = lmer(pl_v_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)
summary(m3_pl_v)
### MASA_d ####
m4_masa = lmer(masa_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)
summary(m4_masa)
### VOL_c ####
m5_vol_c = lmer(vol_c_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)
summary(m5_vol_c)
### log(Area_p) ####
m6_area_p = lmer(area_p_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F) ## log normal ajusta mejor 
summary(m6_area_p)
### l-pappus ####
m7_lpapus = lmer(l_papus_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)  
summary(m7_lpapus)
### l-cipsela ####
m8_lcipsela = lmer(l_cipsela_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)  
summary(m8_lcipsela) ## singular = T , no tiene variación poblacional
### a-cipsela ####
m9_acipsela = lmer(a_cisela_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)  
summary(m9_acipsela) 
### log(ncap) ####
m10_capitulos = lmer(capitulos_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)
summary(m10_capitulos)
### log(altura) ####
m11_altura_planta = lmer(altura_planta_mean ~ 1 +  (1|pob/indiv), data = df_sum_f, REML = F)
summary(m11_altura_planta)

## Tabla y grafico partición de la varianza ####
### tabla ####
table_var_dist = data.frame(
  traits = c('SV', 'PL_m', 'PL_v', 'masa_d', 'vol_c', 'log_area_p', 'l_pappus', 'l_cipsela', 'a_cipsela', 'log_ncap', 'log_alt_p'),
  mean = c(mean(df_sum_f$sv_mean),mean(df_sum_f$pl_m_mean, na.rm = T), mean(df_sum_f$pl_v_mean), mean(df_sum_f$masa_mean, na.rm = T),
           mean(df_sum_f$vol_c_mean), mean(log(df_sum_f$area_p_mean)), mean(df_sum_f$l_papus_mean), mean(df_sum_f$l_cipsela_mean), 
           mean(df_sum_f$a_cisela_mean), mean(log(df_sum_f$capitulos_mean)), mean(log(df_sum_f$altura_planta_mean))),
  sd = c(sd(df_sum_f$sv_mean, na.rm = T), sd(df_sum_f$pl_m_mean, na.rm = T), sd(df_sum_f$pl_v_mean, na.rm = T), 
         sd(df_sum_f$masa_mean, na.rm = T), sd(df_sum_f$vol_c_mean, na.rm = T), sd(log(df_sum_f$area_p_mean), na.rm = T),
         sd(df_sum_f$l_papus_mean, na.rm = T), sd(df_sum_f$l_cipsela_mean, na.rm = T),sd(df_sum_f$a_cisela_mean, na.rm = T),
         sd(log(df_sum_f$capitulos_mean), na.rm = T), sd(log(df_sum_f$altura_planta_mean), na.rm = T)),
  var.pob = c(0.053, 0.0007, 0.0009, 0.0864, 0.1352, 30.81, 0.130, 0, 0.04, 0.7831, 0.1159),
  var.indiv = c(0.1197, 0.0026,0.0016, 0.4003, 0.238, 920.11, 0.652, 0.091, 0.06,1.0772, 0.190248),
  var.res = c(0.1171, 0.0034, 0.0019, 0.49, 0.2869, 228.16, 0.33, 0.05, 0.0912,0.02263, 0.006168)
)
table_var_dist$cv = table_var_dist$sd/table_var_dist$mean

#table_var_dist2 = data.frame(
##  traits = c('l_pappus', 'log_area_p', 'l_cipsela','a_cipsela', 'vol_c','masa_d', 'PL_m', 'PL_v' , 'log_ncap', 'log_alt_p'),
#  mean = c( mean(df_sum_f$l_papus_mean),mean(log(df_sum_f$area_p_mean)),mean(df_sum_f$l_cipsela_mean), mean(df_sum_f$a_cisela_mean),mean(df_sum_f$vol_c_mean), mean(df_sum_f$masa_mean, na.rm = T), 
#            mean(df_sum_f$pl_m_mean, na.rm = T), mean(df_sum_f$pl_v_mean), mean(log(df_sum_f$capitulos_mean)), mean(log(df_sum_f$altura_planta_mean))),
#  sd = c(sd(df_sum_f$l_papus_mean, na.rm = T), sd(log(df_sum_f$area_p_mean), na.rm = T), sd(df_sum_f$l_cipsela_mean, na.rm = T), sd(df_sum_f$a_cisela_mean, na.rm = T), sd(df_sum_f$vol_c_mean, na.rm = T),  sd(df_sum_f$masa_mean, na.rm = T),
#         sd(df_sum_f$pl_m_mean, na.rm = T), sd(df_sum_f$pl_v_mean, na.rm = T), 
#         sd(log(df_sum_f$capitulos_mean), na.rm = T), sd(log(df_sum_f$altura_planta_mean), na.rm = T)),
#  var.pob = c(0.130, 920.11,0,0.04,0.1352, 0.0864, 0.0007, 0.0009, 0.7831, 0.1159),
#  var.indiv = c(0.652, 30.81, 0.091,0.06,0.238,0.4003,0.0026,0.0016, 1.0772, 0.190248),
#  var.res = c(0.33, 228.16, 0.05, 0.0912, 0.2869,0.49,0.0034, 0.0019, 0.02263, 0.006168)
#)

table_var_dist2 = data.frame(
  traits = c('area p','vol c','masa d', 'n cap', 'altura p'),
  mean = c(mean(df_sum_f$area_p_mean), mean(df_sum_f$vol_c_mean), mean(df_sum_f$masa_mean, na.rm = T), 
            mean(df_sum_f$capitulos_mean), mean(df_sum_f$altura_planta_mean)),
  sd = c(sd(df_sum_f$area_p_mean, na.rm = T), sd(df_sum_f$vol_c_mean, na.rm = T), sd(df_sum_f$masa_mean, na.rm = T),
         sd(df_sum_f$capitulos_mean, na.rm = T), sd(df_sum_f$altura_planta_mean, na.rm = T)),
  var.pob = c(30.81, 0.01828, 0.007469, 46233.50, 87.2790),
  var.indiv = c(920.11, 0.05666, 0.160239, 376031.36, 261.7195),
  var.res = c(228.16, 0.08229,0.240154, 23.87, 0.4244)
)

table_var_dist2$cv = table_var_dist2$sd/table_var_dist2$mean

### preparación tabla ####
var_dist_melt = melt(table_var_dist[,c(1,4:6)], id = c('traits')) 
var_dist_melt2 = melt(table_var_dist2[,c(1,4:6)], id = c('traits')) 

## calculo el porentage para cada varianza dentro de cada triat
var_dist_melt <- var_dist_melt %>%  
  group_by(traits) %>%
  mutate(percentage = value / sum(value) * 100)
var_dist_melt$traits = factor(var_dist_melt$traits, levels = (c('log_alt_p','log_ncap', 'l_pappus', 'log_area_p', 'l_cipsela', 'a_cipsela', 'vol_c', 'masa_d','PL_v', 'PL_m','SV'))) ## ordeno
var_dist_melt$variable = factor(var_dist_melt$variable, labels = c('Popultation', 'Individual', 'Residual'))

var_dist_melt2 <- var_dist_melt2 %>%  
  group_by(traits) %>%
  mutate(percentage = value / sum(value) * 100)
var_dist_melt2$traits = factor(var_dist_melt2$traits, levels = (c('altura p', 'n cap', 'masa d', 'vol c', 'area p'))) ## ordeno
var_dist_melt2$traits = factor(var_dist_melt2$traits, labels = (c('altura planta', 'n capítulos', 'masa diáspora', 'vol cipsela', 'área papus')))
var_dist_melt2$variable = factor(var_dist_melt2$variable, labels = c('Localidad', 'Individual', 'Residual + intra-individuo'))
### gráfico partición de la varianza (%) ####
ggplot(var_dist_melt, aes(fill = variable, y = percentage, x = traits)) + 
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = sprintf("%.0f%%", percentage)), position = position_fill(vjust = 0.5)) +
  scale_fill_viridis(discrete = TRUE) + coord_flip() +
  labs(y = "Variance Percentage", x = "Traits")  + my_theme

var_dist_melt2$percentage[13] = 58
ggplot(var_dist_melt2, aes(fill = variable, y = percentage, x = traits)) + 
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = sprintf("%.0f%%", percentage)), position = position_fill(vjust = 0.5)) +
  scale_fill_viridis(discrete = TRUE) + coord_flip() +
  labs(y = "Porcentaje de varianza explicada", x = "Rasgos")  + my_theme
### Grafico partición solo para SV
var_SV = var_dist_melt %>% filter(traits == 'SV')
v_plot = ggplot(var_SV, aes(fill = variable, y = percentage, x = traits)) + 
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = sprintf("%.0f%%", percentage)), position = position_fill(vjust = 0.5)) +
  scale_fill_viridis(discrete = TRUE) + scale_y_continuous(expand = expansion()) + 
  scale_x_discrete(expand = expansion()) + coord_flip() +
  labs(y = "Variance Percentage", x = "Traits")  + my_theme + theme(axis.text.x = element_blank(),
                                                                    axis.ticks = element_blank(),
                                                                    axis.title.y = element_blank(),
                                                                    panel.grid.major.y = element_blank())

plots = wrap_plots(p, v_plot, nrow = 2, heights = c(4,1))
p + plots
plots = wrap_plots(p + ggtitle("a)"), p2 + ggtitle("b)"), ncol = 2, widths = c(1,1))
wrap_plots(plots, v_plot + ggtitle("c)"), nrow = 2, heights = c(5,1))

## graficos de sv en funcion de variables de morfología del fruto ####
### PL_m ####
p_plmsv <- ggplot(df_sum, aes(x = pl_m_mean, y = sv_mean))
p_plmsv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### PL_v ####
# feisimo
p_plvsv <- ggplot(df_sum, aes(x = pl_v_mean, y = sv_mean))
p_plvsv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### masa ####
p_msv <- ggplot(df_sum, aes(x = masa_mean, y = sv_mean))
p_msv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### vol ####
# muy feo
p_vsv <- ggplot(df_sum, aes(x = vol_c_mean, y = sv_mean))
p_vsv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme

## graficos agrupando a nivel de individuo ####
### tabla ####
df_ind <- df2 %>% dplyr::select(pob, indiv, sv, pl_m, pl_v, masa, l_cipsela, a_cisela, l_papus, area_p, vol_c, capitulos, altura_planta ) %>% ## selecciono varaibles de interes
  group_by(pob,indiv)  %>% summarise_all(lst(mean, sd))# marcamos las variables que agrupan
df_ind
### distribución sv ####
p <- ggplot(df_ind, aes(x = sv_mean))
p <- p + geom_histogram(aes(x = sv_mean, y = after_stat(density)),
                        col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                        fill = 'white', bins = 30) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) +
  xlab('Velocidad de Asentamiento media por indiv (m/s)') + ylab('Densidad') +
  scale_y_continuous(expand = c(0,0))+ my_theme
p
p2 <- ggplot(df_ind, aes(x = sv_sd)) ### error
p2 <- p2 + geom_histogram(aes(x = sv_sd, y = ..density..),
                          col = viridis(1, begin = 0.3, option = 'A', direction = 1),
                          fill = 'white', bins = 30) + 
  geom_density(color = viridis(1, begin = 0.3, option = 'A', direction = 1), 
               lwd = 1, 
               fill = viridis(1, begin = 0.3, option = 'A', direction = 1),
               alpha = 0.3) +
  xlab('sd de Velocidad de Asentamiento (m/s)') + ylab('Densidad') +
  scale_y_continuous(expand = c(0,0)) + my_theme
p2
mean(print(df_ind$sv_sd/df_ind$sv_mean))  ### logicamente es mas alta pero no es tan alta
### sv~variables morfologicas ####
### PL_m ####
p_plmsv <- ggplot(df_ind, aes(x = pl_m_mean, y = sv_mean))
p_plmsv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### PL_v ####
# feisimo
p_plvsv <- ggplot(df_ind, aes(x = pl_v_mean, y = sv_mean))
p_plvsv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### masa ####
p_msv <- ggplot(df_ind, aes(x = masa_mean, y = sv_mean))
p_msv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### vol ####
# muy feo
p_vsv <- ggplot(df_ind, aes(x = vol_c_mean, y = sv_mean))
p_vsv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### capitulos ####
# muy feo
head(df_ind)
p_csv <- ggplot(df_ind, aes(x = capitulos_mean, y = sv_mean))
p_csv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### altura de la planta ####
p_asv <- ggplot(df_ind, aes(x = altura_planta_mean, y = sv_mean))
p_asv + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme

### variables morfologicas entre si ####
### masa~areap ####
p_map <- ggplot(df_ind, aes(x = masa_mean, y = area_p_mean))
p_map + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### vol ~ areap
p_vap <- ggplot(df_ind, aes(x = vol_c_mean, y = area_p_mean))
p_vap + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### masa ~ ncap
p_mnc <- ggplot(df_ind, aes(x = capitulos_mean, y = masa_mean))
p_mnc + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
### pl_m ~ ncap
p_plmc <- ggplot(df_ind, aes(x = capitulos_mean, y = pl_m_mean))
p_plmc + geom_smooth(method = 'lm',aes(color = pob), fill = "transparent") +
  geom_smooth(method = 'lm',color = 'black', fill = "transparent", linetype = 'dashed') +  
  geom_point(aes(color = pob)) + scale_color_viridis(discrete = TRUE) + my_theme
## Independent effects analysis ####
sum_df_ord = df_sum_f[, c(5,7,6,10,11,8,12,9)]
hier.part(df_sum_f$sv_mean, sum_df_ord, fam = "gaussian", link = "identity", gof = 'Rsqu')
## modelo sv ~ variables morfologicas ####
model_f = lmer(sv_mean ~ scale(pl_m_mean) + scale(pl_v_mean) + scale(l_papus_mean) + scale(capitulos_mean) + scale(altura_planta_mean) + (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
summary(model_f) ### mucha diferencia en maginitud de efecto entre pl_m vs pl_v
anova(model_f) 
## significancia de las variables  ####
### pl_m ####
model_f1 =lmer(sv_mean ~ scale(pl_v_mean) + scale(l_papus_mean) + scale(capitulos_mean) + scale(altura_planta_mean) + (1|pob/indiv), data = na.omit(df_sum_f), REML = T) 
anova(model_f1,model_f) ## obviamente significativa
### pl_v ####
model_f2 = lmer(sv_mean ~ scale(pl_m_mean) + scale(l_papus_mean) + scale(capitulos_mean) + scale(altura_planta_mean) + (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
anova(model_f2,model_f) ## significativa
### l_papus ####
model_f3 = lmer(sv_mean ~ scale(pl_m_mean) + scale(pl_v_mean) + scale(capitulos_mean) + scale(altura_planta_mean) + (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
anova(model_f3,model_f) ### no significativa
### capitulos ####
model_f4 = lmer(sv_mean ~ scale(pl_m_mean) + scale(pl_v_mean) + scale(l_papus_mean) + scale(altura_planta_mean) + (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
anova(model_f4,model_f) ### no significativa
### alt_planta ####
model_f5 = lmer(sv_mean ~ scale(pl_m_mean) + scale(pl_v_mean) + scale(l_papus_mean) + scale(capitulos_mean) + (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
anova(model_f5,model_f) ### no siginificativo

### modelo pappus area vs masa 
head(df_sum_f)
model_masa = lmer(area_p_mean ~ masa_mean + (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
summary(model_masa)
model_masa_2 = lmer(area_p_mean ~  (1|pob/indiv), data = na.omit(df_sum_f), REML = T)
anova(model_masa_2,model_masa) ## no es significativa la masa para explicar el area del papus, esto quiere decir que
## un aumento de masa no esta compensado con un aumento en el area del pappus por ende seran peores dispersados,
## fortalece la idea de un trade off entre dispersión y dormancia

## distribucion de la varianza juntando los cuadrados ####
### tabla ####
cb <- df2 %>% filter(pob == 'cb') 
cm <- df2 %>% filter(pob == 'cm')
ca <- df2 %>% filter(pob == 'ca')
resto <- df2 %>% filter(pob != 'cb' & pob != 'cm' & pob != 'ca')
unique(cb$indiv)
unique(cm$indiv)
unique(ca$indiv)
cm$indiv[cm$indiv == 1] = 12
cm$indiv[cm$indiv == 2] = 13
cm$indiv[cm$indiv == 3] = 14
cm$indiv[cm$indiv == 4] = 15
cm$indiv[cm$indiv == 5] = 16
cm$indiv[cm$indiv == 6] = 17
cm$indiv[cm$indiv == 7] = 18 
cm$indiv[cm$indiv == 9] = 19
cm$indiv[cm$indiv == 10] = 20
ca$indiv[ca$indiv == 1] = 21
ca$indiv[ca$indiv == 2] = 22
ca$indiv[ca$indiv == 3] = 23
ca$indiv[ca$indiv == 4] = 24
ca$indiv[ca$indiv == 5] = 25
ca$indiv[ca$indiv == 6] = 26
ca$indiv[ca$indiv == 7] = 27 
ca$indiv[ca$indiv == 8] = 28
ca$indiv[ca$indiv == 9] = 29
ca$indiv[ca$indiv == 10] = 30

df3 <- rbind(cb,cm,ca,resto)


for (i in 1:length(df3[,'pob'])) {
  if (df3[i,'pob'] == 'cb' || df3[i,'pob'] == 'cm' || df3[i,'pob'] == 'ca') {
    df3[i,'pob'] <- 'cuadrado'
  }
}

df_sum_f2 <- df3 %>% dplyr::select(pob, indiv, fruto, sv, pl_m, pl_v, masa, l_cipsela, a_cisela, l_papus, area_p, vol_c, capitulos, altura_planta ) %>% ## selecciono varaibles de interes
  group_by(pob,indiv,fruto)  %>% summarise_all(lst(mean))# marcamos las variables que agrupan
df_sum_f2
df_sum_f2$pob = as.factor(df_sum_f2$pob)
df_sum_f2$indiv = as.factor(df_sum_f2$indiv)
### modelos ####
### SV ####
m1_sv2 = lmer(sv_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = T)
summary(m1_sv2) ## la residual incluyo intra y entre frutos
### PL_m ####
m2_pl_m2 = lmer(pl_m_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)
summary(m2_pl_m2)
### PL_v ####
m3_pl_v2 = lmer(pl_v_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)
summary(m3_pl_v2)
### MASA_d ####
m4_masa2 = lmer(masa_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)
summary(m4_masa2)
### VOL_c ####
m5_vol_c2 = lmer(vol_c_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)
summary(m5_vol_c2)
### log(Area_p) ####
m6_area_p2 = lmer(area_p_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = T) ## log normal ajusta mejor 
summary(m6_area_p2)
### l-pappus ####
m7_lpapus2 = lmer(l_papus_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)  
summary(m7_lpapus2)
### l-cipsela ####
m8_lcipsela2 = lmer(l_cipsela_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)  
summary(m8_lcipsela2) ## singular = T , no tiene variación poblacional
### a-cipsela ####
m9_acipsela2 = lmer(a_cisela_mean ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)  
summary(m9_acipsela2) 
### log(ncap) ####
m10_capitulos2 = lmer(log(capitulos_mean) ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)
summary(m10_capitulos2)
### log(altura) ####
m11_altura_planta2 = lmer(log(altura_planta_mean) ~ 1 +  (1|pob/indiv), data = df_sum_f2, REML = F)
summary(m11_altura_planta2)

