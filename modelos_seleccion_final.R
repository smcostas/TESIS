library(tidyverse);library(viridis);library(lme4); library(fitdistrplus); library(MuMIn); library(cAIC4); library(parallel);library(pbkrtest)
library(ggcorrplot); library(PerformanceAnalytics); library(reshape2); library(hier.part); library(patchwork); library(ggeffects)
library(mgcv);library(boot); library(parallel); library(mgcv);library(gamm4)
library(gratia);library(tidygam); library(metR)

# ARMADO DE TABLAS ####################
sub_df = read.csv("sub_df2.csv") ### esto si quiero usar la clz obtenida a partir de la estandarizacion del settling velocity
#sub_df$clz = sub_df$inverse_sv*sub_df$prob_germ #esto si quiero usar la clz obtenida con sv sin estandarizar
sub_df <- sub_df %>%
  group_by(pob, indiv) %>% # Agrupamos por población e individuo
  mutate(indiv2 = cur_group_id()) %>% # Asignamos un ID único continuo
  ungroup()
sub_df_fr <- sub_df %>%
  dplyr::select(pob, indiv, indiv2,fruto, s_inverse_sv, prob_germ,germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta) %>%
  group_by(pob, indiv,fruto) %>%
  summarise(across(c(indiv2, s_inverse_sv, prob_germ,germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta), mean, .names = "{col}"))

mean(sub_df_fr$masa)
sub_df_fr$W = sub_df_fr$clz/mean(sub_df_fr$clz)
max(sub_df_fr$clz)
summary(sub_df_fr$W)
summary(sub_df_fr$s_inverse_sv)
#plot(sub_df_fr$clz ~ sub_df_fr$masa)
#plot(sub_df_fr$W ~ sub_df_fr$masa)
corr <- round(cor(na.omit(sub_df_fr[,c('masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')])), 3)

corr
corrplot = ggcorrplot(corr, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T, legend.title = 'r')

corrplot = corrplot + scale_x_discrete(labels = c('area_p' = 'área papus',
                                       'vol_c' = 'volumen Cipsela',
                                       'capitulos' = 'capítulos',
                                       'altura_planta' = 'altura planta')) +
            scale_y_discrete(labels = c('masa' = 'masa diáspora',
                              'area_p' = 'área papus',
                              'vol_c' = 'volumen cipsela',
                              'capitulos' = 'capítulos'
                              ))
### grafico corrplot#######
pdf('corrplot.pdf')
corrplot
dev.off()
png('corrplot.png',width = 15, height = 10, units = 'cm', res = 600)
corrplot  
dev.off()

# SELECCION FENOTIPICA A NIVEL FILIAL ############
### estandarizacion #######
est_df = sub_df_fr
est_v = c('masa', 'area_p', 'vol_c', 'capitulos', 
                        'altura_planta', 'pl_m', 'pl_v')

min(sub_df_fr$masa)

for (col in est_v) {
  est_df[,col] <- as.vector(scale(est_df[,col]))
}

## OPORTUNIDAD DE SELECCION ###############
w.os = var(sub_df_fr$W) ## var 0.4349
sd(sub_df_fr$W)/mean(sub_df_fr$W) ### coeficiente de variacion alto 0.659
w.os/sqrt(2*length(sub_df_fr$W)) # error estandar de 0.0232

## DIFERENCIAL DE SELECCION LINEAL #######################
fixef_fun <- function(model) {
  return(fixef(model))
}
n_cores <- detectCores() - 1
conf_levels <- c(0.90, 0.95, 0.99, 0.999)

### masa ####
ds_masa = lmer(W ~ masa + (1|indiv2), data = est_df)
summary(ds_masa)
boot_results <- bootMer(ds_masa, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)

ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##***
### area_p #############
ds_areap = lmer(W ~ area_p + (1|indiv2), data = est_df)
summary(ds_areap)
boot_results <- bootMer(ds_areap, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##***

### vol #########
ds_volc = lmer(W ~ vol_c + (1|indiv2), data = est_df)
ds_volc2 = lmer(W ~ (1|indiv2), data = est_df)

boot_results <- bootMer(ds_volc, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##**


## diferencial de seleccion no lineal ########
## masa ########
ds_masa2 = lmer(W ~ I(0.5*(masa^2)) + (1|indiv2), data = est_df)
summary(ds_masa2)

boot_results <- bootMer(ds_masa2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## significativo
print(ci)
## areap #######
ds_areap2 = lmer(W ~ I(0.5*(area_p^2)) + (1|indiv2), data = est_df)
summary(ds_areap2)

boot_results <- bootMer(ds_areap2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### vol #########
ds_volc2 = lmer(W ~  I(0.5*(vol_c^2)) + (1|pob/indiv), data = est_df)
summary(ds_volc2)
boot_results <- bootMer(ds_volc2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

beep()
## GRADIENTE DE SELECCION DIRECCIONAL ##########################
### modelo ###################
grs = lmer(W ~ area_p + masa + vol_c + (1|indiv2), data = est_df)
summary(grs)

### bootstrap para significancia
boot_results <- bootMer(grs, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)

### intervalos de confianza #######

get_conf_intervals <- function(results, conf_levels) {
  conf_intervals <- list()
  for (i in 1:ncol(results$t)) {
    ci <- boot.ci(results, index = i, type = "perc", conf = conf_levels)
    conf_intervals[[names(results$t0)[i]]] <- ci ## nombro las sub listas como la variable
  }
  return(conf_intervals)
}
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(boot_results, conf_levels)


## GRADIENTE de SELECCION CUADRATICA y CORRELACIONAL ###########
### modelo #########################
grs_cc = lmer(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
                I(0.5*(vol_c^2)) + area_p:masa + masa:vol_c + area_p:vol_c  + 
                (1|indiv2), data = est_df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(grs_cc)
#grs_cc = lm(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
#                I(0.5*(vol_c^2)) + area_p:masa + masa:vol_c + area_p:vol_c, data = est_df)

## creo la funcion necesaria para para remuestrear y obtener la distribucon del modelo

## otra forma es coon butmer de lme4 que es especifico para modelos mixtos, respeta la estructura del moel ## 

### si no usara efecto aleatroio
#coef_fun <- function(data, indices) {
#
#  d <- data[indices, ]
#
#  model <- lmer(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) + 
#                  I(0.5*(area_p^2)) + I(0.5*(vol_c^2)) + 
#                  area_p:masa + masa:vol_c + area_p:vol_c + 
#                  (1|pob/indiv), data = d, 
#                control = lmerControl(optimizer = "bobyqa", 
#                                      optCtrl = list(maxfun = 100000)))
#  return(fixef(model))  ## obtengo los coeficientes del modelo en cada corrida
#}
#coef_fun <- function(data, indices) {
#  
#  d <- data[indices, ]
#  
#  model <- lm(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) + 
#                I(0.5*(area_p^2)) + I(0.5*(vol_c^2)) + 
#                area_p:masa + masa:vol_c + area_p:vol_c, data = d)
#  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
#}

### boostrap para significancia ###################
boot_results <- bootMer(grs_cc, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
#n_cores <- detectCores() - 1
#cl <- makeCluster(n_cores)
#set.seed(123) 
#results <- boot(data = est_df, statistic = coef_fun, R = 10000)
#stopCluster(cl)

### intervalos de confianza #######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(boot_results, conf_levels)




library(plotly)
##  efectos de la seleccion cuadratica y correlacional ########

grid_data <- expand.grid(
  masa = seq(min(est_df$masa), max(est_df$masa), length = 100),
  area_p = seq(min(est_df$area_p), max(est_df$area_p), length = 100)
)
grid_data$vol_c = 0 ## para el volumen en la media
grid_data$predicted <- predict(grs_cc, newdata = grid_data, re.form = NA)


p <- ggplot(grid_data, aes(x = masa, y = area_p, z = predicted)) +
  geom_tile(aes(fill = predicted)) +
  geom_contour(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Masa de la diáspora", y = "Área del papus", fill = "Colonización")
p
plotly::ggplotly(p)
masa_seq = seq(min(est_df$masa), max(est_df$masa), length = 100)
area_p_seq = seq(min(est_df$area_p), max(est_df$area_p), length = 100)
predicted <- matrix(grid_data$predicted, nrow = length(masa_seq), ncol = length(area_p_seq))
plot_ly(
  x = ~masa_seq,
  y = ~area_p_seq,
  z = ~predicted,
  type = "surface",
  colorscale = 'Viridis',
  colorbar = list(title = "Colonización\npredicha")
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Masa"),
      yaxis = list(title = "Área del papus"),
      zaxis = list(title = "Colonización predicha")
    )
  )



## volumen y area
grid_data <- expand.grid(
  vol_c = seq(min(est_df$vol_c), max(est_df$vol_c), length = 100),
  area_p = seq(min(est_df$area_p), max(est_df$area_p), length = 100)
)
grid_data$masa = 0 ## para el volumen en la media
grid_data$predicted <- predict(grs_cc, newdata = grid_data, re.form = NA)


p <- ggplot(grid_data, aes(x = vol_c, y = area_p, z = predicted)) +
  geom_tile(aes(fill = predicted)) +
  geom_contour(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "gradiente de Selección correlacional y no lineal", x = "Vol cipsela", y = "Área de papus", fill = "Fitness")
p
plotly::ggplotly(p)


plot_ly(
  x = ~grid_data$vol_c,
  y = ~grid_data$area_p,
  z = ~grid_data$predicted,
  type = "scatter3d",
  mode = "lines",
  line = list(width = 3, color = ~grid_data$predicted, colorscale = 'Viridis')
) %>%
  layout(
    title = "Gradiente de seleccion correlacional y no lineal",
    scene = list(
      xaxis = list(title = "Volumen cipsela"),
      yaxis = list(title = "area del pappus"),
      zaxis = list(title = "Fitness")
    )
  )



## GAMs gradientes no lineales ####
plot_model = function(model, x, y, xname) {
  newx = seq(min(x), max(x), length.out = 100)
  
  newdata = data.frame(newx)
  colnames(newdata) <- xname
  
  newy = predict(model, newdata = newdata, se.fit = TRUE, exclude = "s(indiv2)",newdata.guaranteed=TRUE, type = 'response')
  
  yhat = newy$fit
  yup = newy$fit + newy$se.fit
  ydown = newy$fit - newy$se.fit
  
  df = data.frame(yhat = yhat, newx = newx, yup = yup, ydown = ydown)
  
  p <- ggplot(df, aes(x = newx, y =  yhat)) 
  p <- p + geom_ribbon(aes(ymin = ydown, ymax = yup), 
                        alpha = .15, color = 'grey', fill = 'grey', linetype = 'dashed') + geom_line(aes(y = yhat), 
                        linewidth = 0.3) + xlab(xname)  + ylab('colonización')
  return(p)
}
### masa ########
sub_df_fr$indiv2 <- as.factor(sub_df_fr$indiv2)
fit = gam(W ~ s(masa, bs = 'cr') + s(indiv2, bs = 're') , data = sub_df_fr, method = 'REML')
p = plot_model(fit, sub_df_fr$masa, sub_df_fr$W,'masa')
p_m = p + geom_point(data = sub_df_fr, aes(x = masa , y = W), color = 'grey14', alpha = 0.3) + 
  scale_x_continuous(breaks = seq(0,3,0.5)) + xlab('masa (mg)') + theme_bw()
p_m = p_m + geom_vline(aes(xintercept= mean(sub_df_fr$masa)), 
                 linetype = 'dashed', color = 'green') 
### papus ###### 
fit = gam(W ~ s(area_p, bs = 'cr')+ s(indiv2, bs = 're') , data = sub_df_fr, method = 'REML')
p = plot_model(fit, sub_df_fr$area_p, sub_df_fr$W,'area_p')
p_p = p + geom_point(data = sub_df_fr, aes(x = area_p , y = W, color = prob_germ), alpha = 0.7) + 
  labs(x = 'área del papus (mm2)', color = 'prob\ngerminación') + 
  scale_x_continuous(breaks = seq(100,250,25)) + scale_color_viridis_c(option = 'A') + theme_bw()
p_p = p_p + geom_vline(aes(xintercept= mean(sub_df_fr$area_p)), 
                 linetype = 'dashed', color = 'green')

eff = wrap_plots(p_m + ggtitle("(a)") +
             p_p + ggtitle("(b)"))


### volumen ######
fit = gam(W ~ s(vol_c, bs = 'cr') + s(indiv2, bs = 're'), data = est_df, method = 'REML')
p = plot_model(fit, est_df$vol_c, est_df$W,'vol_c')
p + geom_point(data = est_df, aes(x = vol_c , y = W), color = 'grey14', alpha = 0.3) + xlab('volumen de la cipsela') + theme_bw()
### graficos final #####
pdf('splines_cubicos_diasporas.pdf')
eff
dev.off()
png('splines_cubicos_diasporas.png',width = 20, height = 15, units = 'cm', res = 600)
eff
dev.off()
## SUPERFICIE DE SELECCION GAMS ####

est_df$indiv2 <- as.factor(est_df$indiv2)
sub_df_fr$indiv2 <- as.factor(sub_df_fr$indiv2)

m3 <- gam(W ~ s(area_p, bs = "cr") + s(masa, bs = "cr" ) + s(indiv2, bs = "re"), 
          data=est_df, method = "REML", 
          select = TRUE)
m4 <- gam(W ~ te(area_p, masa) + s(indiv2, bs = "re"), 
                data = sub_df_fr, method = "REML",
                select = TRUE)
m5 <- gam(W ~ s(area_p, bs = "cr") + s(masa, bs = "cr") + ti(area_p, masa) + s(indiv2, bs = "re"),
          data = est_df, method = "REML", 
          select = TRUE)
m6 = gam(W ~ s(area_p, bs = "cr") + s(masa, bs = "cr") + ti(area_p, masa) + s(indiv2, bs = "re"),
         data = est_df, method = "REML", 
         select = TRUE)
m2 <- gam(W ~ s(area_p, masa, k = 10, bs = "tp") + s(indiv2, bs = "re") , data=sub_df_fr, 
          method = "REML", 
          select = TRUE)
AIC(m3, m4, m5)
summary(m4)
vis.gam(m2, view = c("masa", "area_p"), type = "response", 
        plot.type = "contour", color = "cm")
## esto si quiero usar variables no estandarizadas cuando las use estandarizadas
mu_area <- mean(sub_df_fr$area_p)
sigma_area <- sd(sub_df_fr$area_p)

mu_masa <- mean(sub_df_fr$masa)
sigma_masa <- sd(sub_df_fr$masa)

area_seq <- seq(from = min(sub_df_fr$area_p), to = max(sub_df_fr$area_p), length.out = 50)  # Valores no estandarizados
masa_seq <- seq(from = min(sub_df_fr$masa), to = max(sub_df_fr$masa), length.out = 50)

area_std <- (area_seq - mu_area) / sigma_area
masa_std <- (masa_seq - mu_masa) / sigma_masa
new_data <- expand.grid(area_p = area_std, masa = masa_std)


predictions <- predict(m4, newdata = new_data, exclude = "s(indiv2)",newdata.guaranteed=TRUE)

new_data$area_original <- area_seq[as.numeric(as.factor(new_data$area_p))]
new_data$masa_original <- masa_seq[as.numeric(as.factor(new_data$masa))]
new_data$W_pred <- predictions
## Esto si no estan estandarizadas o si quiero usarlas estandarizadas
new_data = predict_gam(m4, exclude_terms = 's(indiv2)' ,length_out = 50)

ggplot(new_data, aes(x = area_p, y = masa, z = W)) +
  geom_tile(aes(fill = W)) +  # Superficie continua
  geom_contour(aes(z = W),color = 'black', size = 0.3) +  # Contornos blancos
  scale_fill_viridis_c(option = "viridis", name = "Colonización relativa",
                       breaks = seq(0, 2.4, by = 0.2),
                       labels = seq(0, 2.4, by = 0.2)) +  # Escala continua
  labs(x = "Área del Papus (mm2)", 
       y = "Masa de la diáspora (mg)") + scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  theme_bw() + theme(legend.key.height = unit(1.5, "cm"),  # Altura de la barra
                     legend.key.width = unit(0.7, "cm"))
superficie = ggplot(new_data, aes(x = masa, y = area_p, z = W))  + 
  geom_contour_fill(aes(fill = after_stat(level))) + 
  geom_point(data = sub_df_fr, aes(x = masa, y = area_p, z = W), color = 'grey15', alpha = 0.5) +
  scale_fill_discretised(low = "#440154CC", high = "#FDE725CC") + 
  geom_vline(aes(xintercept= mean(sub_df_fr$masa)), 
             linetype = 'dashed', color = 'green') +
  geom_hline(aes(yintercept= mean(sub_df_fr$area_p)), 
             linetype = 'dashed', color = 'green') +
  labs(x = "Masa de la diáspora (mg)", 
       y = "Área del Papus (mm2)",
       fill = "Colonización\n  predicha") + 
  scale_y_continuous(breaks = seq(125, 250, by = 25),
                     limits = c(min(sub_df_fr$area_p), max(sub_df_fr$area_p)),
                     expand = c(0,0)) + 
  scale_x_continuous(breaks = seq(0.5, 2.5, by = .5),
                     limits = c(min(sub_df_fr$masa), max(sub_df_fr$masa)),
                     expand = c(0,0)) +
  theme_bw() + theme(
    legend.key.height = unit(1, "cm"),  # Altura de la barra
    legend.key.width = unit(0.8, "cm"),     # Ancho de la barra
    legend.box = "vertical",            # Poner las leyendas en fila
    legend.position = "right"            # Ubicar las leyendas en la parte inferior
  )
### graficos final #####
pdf('superficie_diasporas.pdf')
superficie  
dev.off()
png('superficie_diasporas.png',width = 20, height = 15, units = 'cm', res = 600)
superficie  
dev.off()

# SELECCION FENOTIPICA NIVEL MATERNO CON TODAS LAS DIASPORAS  ###############################
df_sum <- sub_df %>% dplyr::select(pob, indiv, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta ) %>% group_by(pob,indiv)  %>% summarise_all(lst(mean,sd))# marcamos las variables que agrupan
colnames(df_sum)
sub_df_sum = df_sum[,c(1:12)]

colnames(sub_df_sum) = c('pob', 'indiv', 'sv', 'prob_germ', 'clz', 'pl_m', 'pl_v',
                         'masa diáspora', 'área papus', 'volumen cipsela', 'capítulos', 'altura planta')


r <- round(cor(na.omit(sub_df_sum[,c('masa diáspora', 'área papus', 'volumen cipsela', 'capítulos', 'altura planta')])), 3)


r = ggcorrplot(r, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T, legend.title = 'r') 

pdf('corrplot_indiv.pdf')
r
dev.off()
png('corrplot_indiv.png',width = 15, height = 10, units = 'cm', res = 600)
r  
dev.off()

colnames(sub_df_sum) = c('pob', 'indiv', 'sv', 'prob_germ', 'clz', 'pl_m', 'pl_v',
                         'masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')

est_df = sub_df_sum ## lo nombro igual para facilitar la repeticion
est_v = c('masa', 'area_p', 'vol_c', 'capitulos', 
          'altura_planta', 'pl_m', 'pl_v')

for (col in est_v) {
  est_df[,col] <- as.vector(scale(est_df[,col]))
}

## estandarizacion del fitness relativo. #####################################
sub_df_sum$W  = est_df$W = sub_df_sum$clz/mean(sub_df_sum$clz)
## OPORTUNIDAD DE SELECCION ###############
w.os = var(sub_df_sum$W) ## var 0.269, obviamente disminuye la varianza
sd(sub_df_sum$W)/mean(sub_df_sum$W) ### coeficiente de variacion alto 0.518
w.os/sqrt(2*length(sub_df_sum$W)) # error estandar de 0.0248

## DIFERENCIAL DE SELECCION LINEAL #######################
### masa ####
ds_masa = lmer(W ~ masa + (1|pob), data = est_df)
summary(ds_masa)
ds_1 = lmer(W ~ 1 + (1|pob), data = est_df)
anova(ds_masa, ds_1) ## ahora uso esto para ascelerar tiempos ##  < 2.2e-16 *** masa re contra significativa
### area_p #############
ds_areap = lmer(W ~ area_p + (1|pob), data = est_df)
anova(ds_areap, ds_1) ## significativo 
### vol #########
ds_volc = lmer(W ~ vol_c + (1|pob), data = est_df)
summary(ds_volc)

anova(ds_volc, ds_1) ###### significativo
### capitulos ##################
ds_cap = lmer(W ~ scale(capitulos) + (1|pob), data = est_df)
summary(ds_cap)
anova(ds_cap, ds_1) ## no significativo
### altura planta ######
ds_altp = lmer(W ~ altura_planta + (1|pob), data = est_df)
summary(ds_altp)
anova(ds_altp, ds_1)  ## signficativa (marginalmente 0.04949):O


## diferencial de seleccion no lineal ########
### masa ########
ds_masa2 = lmer(W ~ I(0.5*(masa^2)) + (1|pob), data = est_df)
summary(ds_masa2)
anova(ds_masa2, ds_1) ## el componente cuadratico es no significativo
### areap #######
ds_areap2 = lmer(W ~ I(0.5*(area_p^2)) + (1|pob), data = est_df)
summary(ds_areap2)
anova(ds_areap2, ds_1) ### no significativo
### vol #########
ds_volc2 = lmer(W ~  I(0.5*(vol_c^2)) + (1|pob), data = est_df)
anova(ds_volc2, ds_1) ###### no significativo (0.09)
### capitulos ##################
ds_cap = lmer(W ~  I(0.5*(capitulos^2)) + (1|pob), data = est_df)
anova(ds_cap, ds_1) ## no significativo
### altura planta ######
ds_altp = lmer(W ~ I(0.5*(altura_planta^2)) + (1|pob), data = est_df)
anova(ds_altp, ds_1)  ## no signficativa


## GRADIENTE DE SELECCION DIRECCIONAL ##########################
grs = lmer(W ~ area_p + masa + vol_c + capitulos + altura_planta + (1|pob), data = est_df)
summary(grs)
### pappus ########
grs_s1 = lmer(W ~  masa + vol_c + capitulos + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s1) ### significativo 0.00122

### masa ##########
grs_s2 = lmer(W ~  area_p + vol_c + capitulos + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s2) ### significativo p < 2.2e-16 ***

### vol ##########
grs_s3 = lmer(W ~  area_p + masa + capitulos + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s3) ### no significativo 0.1475

### capitulos ##########
grs_s3 = lmer(W ~  area_p + masa + vol_c + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s3) ### no significativo 0.4108

### altura planta ##########
grs_s3 = lmer(W ~  area_p + masa + vol_c + capitulos + (1|pob), data = est_df)
anova(grs,grs_s3) ### no significativo 0.2279

## GRADIENTE de SELECCION CUADRATICA y CORRELACIONAL ########

grs_cc = lm(W ~ area_p + masa + altura_planta + I(0.5*(masa^2)) + I(0.5*(area_p^2)) + I(0.5*(altura_planta^2)) + 
                area_p:masa  + area_p:altura_planta + masa:altura_planta, data = est_df) ## no tiene efecto de la poblacion


summary(grs_cc)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ area_p + masa + altura_planta + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
                I(0.5*(altura_planta^2)) + area_p:masa + masa:altura_planta + area_p:altura_planta, data = d)
  
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
### boostrap #########
cl <- makeCluster(n_cores)
set.seed(123) 
results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
### intervalo sde confianza ######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(results, conf_levels)

## gam #####
fit = gam(W ~ s(masa, bs = 'cr'), data = sub_df_sum, method = 'REML')

p_masa = plot_model(fit, sub_df_sum$masa, sub_df_sum$W,'masa')
p_masa = p_masa + geom_point(data = sub_df_sum, aes(x = masa , y = W), color = 'grey14', alpha = 0.3) + 
  scale_x_continuous(breaks = seq(0,3,0.5)) + xlab('masa (mg)') + theme_minimal()
p_masa = p_masa + geom_vline(aes(xintercept= mean(sub_df_sum$masa)), 
                       linetype = 'dashed', color = 'green') 

fit2 = gam(W ~ s(area_p, bs = 'cr'), data = sub_df_sum, method = 'REML')

p_area = plot_model(fit2, sub_df_sum$area_p, sub_df_sum$W,'area_p')
p_area = p_area + geom_point(data = sub_df_sum, aes(x = area_p , y = W), color = 'grey14', alpha = 0.3) + 
         xlab('Área del papus (mm2)') + theme_minimal()
p_area = p_area + geom_vline(aes(xintercept= mean(sub_df_sum$area_p)), 
                             linetype = 'dashed', color = 'green') 

p_area + p_masa
## SUPERFICIE DE  selección #####
m_indiv <- gam(W ~ te(area_p, masa) , 
                     data = sub_df_sum, method = "REML",
                     select = TRUE)
#m_indiv <- gam(W ~ s(masa, bs = 'cr') + s(area_p, bs = 'cr'), data = sub_df_sum, method = 'REML')


vis.gam(m_indiv, view = c("masa", "area_p"), type = "response", 
        plot.type = "contour", color = "cm")

new_data = predict_gam(m_indiv ,length_out = 50)
plot_superficie = ggplot(new_data, aes(x = masa, y = area_p, z = W))  + 
  geom_contour_fill(aes(fill = after_stat(level))) + 
  geom_point(data = sub_df_sum, aes(x = masa, y = area_p, z = W), color = 'grey15', alpha = 0.5) +
  scale_fill_discretised(low = "#440154CC", high = "#FDE725CC") +  
  geom_vline(aes(xintercept= mean(sub_df_sum$masa)), 
             linetype = 'dashed', color = 'green') +
  geom_hline(aes(yintercept= mean(sub_df_sum$area_p)), 
             linetype = 'dashed', color = 'green') +
  labs(x = "Masa de la diáspora (mg)", 
       y = "Área del Papus (mm2)",
       fill = "Colonización\n  predicha") + 
  scale_y_continuous(breaks = seq(125, 250, by = 25),
                     limits = c(min(sub_df_sum$area_p), max(sub_df_sum$area_p)),
                     expand = c(0,0)) + 
  scale_x_continuous(breaks = seq(0.5, 2.5, by = .5),
                     limits = c(min(sub_df_sum$masa), max(sub_df_sum$masa)),
                     expand = c(0,0)) +
  theme_bw() + theme(
    legend.key.height = unit(0.8, "cm"),  # Altura de la barra
    legend.key.width = unit(0.5, "cm"),     # Ancho de la barra
    legend.box = "vertical",            # Poner las leyendas en fila
    legend.position = "right"            # Ubicar las leyendas en la parte inferior
  )



# SELECCION FENOTIPICA NIVEL MATERNO CON DIASPORAS GERMINADAS UNICAMENTE ####################
var <- sub_df %>%
  filter(germ == 1) %>%
  dplyr::select(pob, indiv, indiv2, fruto, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta) %>%
  group_by(pob, indiv, fruto) %>%
  summarise(across(c(indiv2, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta), mean, .names = "{col}"))


## varianza explicada ####
m1 = lmer(masa ~ (1|indiv2), data = var)
summary(m1)
m2 = lmer(vol_c ~ (1|pob/indiv), data = var)
summary(m2)
m3 = lmer(area_p ~ (1|pob/indiv), data = var)
summary(m3)

## data set final ######
sub_df_sum <- sub_df %>%
  filter(germ == 1) %>%
  dplyr::select(pob, indiv,fruto, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta) %>%
  group_by(pob, indiv) %>%
  summarise(across(c(sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta), mean, .names = "{col}"),
            n = n_distinct(fruto))

mean(sub_df_sum$clz)
var(sub_df_sum$clz)

corr <- round(cor(na.omit(sub_df_sum[,c('masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')])), 3)

corr
corrplot = ggcorrplot(corr, type = 'lower',outline.color = "white",
                      ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T, legend.title = 'r')

corrplot = corrplot + scale_x_discrete(labels = c('area_p' = 'área papus',
                                                  'vol_c' = 'volumen Cipsela',
                                                  'capitulos' = 'capítulos',
                                                  'altura_planta' = 'altura planta')) +
  scale_y_discrete(labels = c('masa' = 'masa diáspora',
                              'area_p' = 'área papus',
                              'vol_c' = 'volumen cipsela',
                              'capitulos' = 'capítulos'
  ))
### grafico corrplot#######
#pdf('corrplot_indivgerm.pdf')
#corrplot
#dev.off()
#png('corrplot_indivgerm.png',width = 15, height = 10, units = 'cm', res = 600)
#corrplot  
#dev.off()
est_df = sub_df_sum ## lo nombro igual para facilitar la repeticion
est_v = c('masa', 'area_p', 'vol_c', 'capitulos', 
          'altura_planta', 'pl_m', 'pl_v')

for (col in est_v) {
  est_df[,col] <- as.vector(scale(est_df[,col]))
}

## estandarizacion del fitness relativo. #####################################
sub_df_sum$W  = est_df$W = sub_df_sum$clz/mean(sub_df_sum$clz)
## OPORTUNIDAD DE SELECCION ###############
w.os = var(sub_df_sum$W) ## var 0.031, obviamente disminuye la varianza
sd(sub_df_sum$W)/mean(sub_df_sum$W) ### coeficiente de variacion alto 0.17
w.os/sqrt(2*length(sub_df_sum$W)) # error estandar de 0.003065849

## efecto de la poblacion para simplificar #######
grs_cc = lmer(W ~ area_p + masa + vol_c + altura_planta + capitulos  + I(0.5*(masa^2)) + I(0.5*(area_p^2)) + 
              I(0.5*(vol_c^2)) + I(0.5*(capitulos^2)) + I(0.5*(altura_planta^2)) + 
              area_p:masa + masa:vol_c + area_p:vol_c + area_p:capitulos + area_p:altura_planta + 
              masa:capitulos + masa:altura_planta + vol_c:capitulos + 
              vol_c:altura_planta + capitulos:altura_planta + (1|pob),data = est_df)
grs_ccs = lm(W ~ area_p + masa + vol_c + altura_planta + capitulos  + I(0.5*(masa^2)) + I(0.5*(area_p^2)) + 
                I(0.5*(vol_c^2)) + I(0.5*(capitulos^2)) + I(0.5*(altura_planta^2)) + 
                area_p:masa + masa:vol_c + area_p:vol_c + area_p:capitulos + area_p:altura_planta + 
                masa:capitulos + masa:altura_planta + vol_c:capitulos + 
                vol_c:altura_planta + capitulos:altura_planta ,data = est_df)

anova(grs_cc, grs_ccs) ## no hay efecto de la poblacion


  

## DIFERENCIAL DE SELECCION LINEAL #######################
### masa ####
ds_masa = lm(W ~ masa, data = est_df)
#summary(ds_masa)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ masa, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##  significativo
print(ci)
### area_p #############
ds_areap = lm(W ~ area_p, data = est_df)
#summary(ds_areap)## significativo

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ area_p, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## significativo
print(ci)

### vol #########
ds_volc = lm(W ~ vol_c, data = est_df)
#summary(ds_volc)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ vol_c, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### capitulos ##################
ds_cap = lm(W ~ capitulos , data = est_df)
#summary(ds_cap)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ capitulos, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### altura planta ######
ds_altp = lm(W ~ altura_planta, data = est_df)
#summary(ds_altp)  ## 
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ altura_planta, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

beep()
## diferencial de seleccion no lineal ########
### masa ########
ds_masa2 = lm(W ~ I(0.5*(masa^2)), data = est_df)
#summary(ds_masa2)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(masa^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## **significativo
print(ci)
 ## el componente cuadratico es significativo 0.002393
### areap #######
ds_areap2 = lm(W ~ I(0.5*(area_p^2)), data = est_df)
#summary(ds_areap2)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(area_p^2)), data = d)
  return(coef(model))  
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

### no significativo
### vol #########
ds_volc2 = lm(W ~  I(0.5*(vol_c^2)) , data = est_df)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(vol_c^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### capitulos ##################
ds_cap2 = lm(W ~  I(0.5*(capitulos^2)), data = est_df)
#summary(ds_cap2) ## no significativo
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(capitulos^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## *significativo
print(ci)
### altura planta ######
ds_altp2 = lm(W ~ I(0.5*(altura_planta^2)), data = est_df)
#summary(ds_altp2)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(altura_planta^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

beep()
## GRADIENTE DE SELECCION DIRECCIONAL ##########################
grs = lm(W ~ area_p + masa  + altura_planta , data = est_df)
summary(grs)
### boostrap para significancia #######
coef_fun <- function(data, indices) {
    
    d <- data[indices, ]
    
    model <- lm(W ~ area_p + masa + altura_planta, 
                data = d)
    return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)

### intervalos de confianza #######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(results, conf_levels)

## GRADIENTE de SELECCION CUADRATICA y CORRELACIONAL ########


grs_cc = lm(W ~ area_p + masa + altura_planta + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
              I(0.5*(altura_planta^2)) + area_p:masa + masa:altura_planta + area_p:altura_planta, data = est_df)
summary(grs_cc)
### boostrap para significancia ######

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ area_p + masa + altura_planta + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
                I(0.5*(altura_planta^2)) + area_p:masa + masa:altura_planta + area_p:altura_planta, data = d)
  
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}


cl <- makeCluster(n_cores)
set.seed(123) 
results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
### intervalo sde confianza ######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(results, conf_levels)

##  efectos de la seleccion cuadratica y correlacional ########

area_p_seq <- seq(from = min(est_df$area_p), to = max(est_df$area_p), length.out = 100)
masa_seq <- seq(from = min(est_df$masa), to = max(est_df$masa), length.out = 100)

grid_data <- expand.grid(area_p = area_p_seq, masa = masa_seq, altura_planta = 0)

grid_data$predicted <- predict(grs_cc, newdata = grid_data)

p <- ggplot(grid_data, aes(x = masa, y = area_p, z = predicted)) +
  geom_tile(aes(fill = predicted)) +
  geom_contour(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs( x = "Masa", y = "Área del papus", fill = "Colonización predicha")
p
#pplotly::ggplotly(p)

predicted <- matrix(grid_data$predicted, nrow = length(masa_seq), ncol = length(area_p_seq))

plot_ly(
  x = ~masa_seq,
  y = ~area_p_seq,
  z = ~predicted,
  type = "surface",
  colorscale = 'Viridis',
  colorbar = list(title = "Colonización\npredicha")
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Masa"),
      yaxis = list(title = "Área del papus"),
      zaxis = list(title = "Colonización predicha")
    )
  )
### gams #####

plot_model = function(model, x, y, xname) {
  newx = seq(min(x), max(x), length.out = 100)
  
  newdata = data.frame(newx)
  colnames(newdata) <- xname
  
  newy = predict(model, newdata = newdata, se.fit = TRUE, type = 'response')
  
  yhat = newy$fit
  yup = newy$fit + newy$se.fit
  ydown = newy$fit - newy$se.fit
  
  df = data.frame(yhat = yhat, newx = newx, yup = yup, ydown = ydown)
  
  p <- ggplot(df, aes(x = newx, y =  yhat)) 
  p <- p + geom_ribbon(aes(ymin = ydown, ymax = yup), 
                       alpha = .15, color = 'grey', 
                       fill = 'grey', linetype = 'dashed') + 
            geom_line(aes(y = yhat), linewidth = 0.3) + xlab(xname)  + 
    ylab('colonización')
  return(p)
}

### masa ########
fit = gam(W ~ s(masa, bs = 'cr'), data = sub_df_sum, method = 'REML')
p = plot_model(fit, sub_df_sum$masa, sub_df_sum$W,'masa')
p_m = p + geom_point(data = sub_df_sum, aes(x = masa , y = W), 
                     color = 'grey14', alpha = 0.3) + scale_x_continuous(breaks = seq(0,3,0.5)) + 
          xlab('masa (mg)') + theme_bw() 
p_m = p_m + geom_vline(aes(xintercept = mean(sub_df_sum$masa)), 
                 linetype = 'dashed', color = 'green')

### papus ###### 
#sub_df_sum  = sub_df_sum[sub_df_sum$area_p < max(sub_df_sum$area_p),]
fit = gam(W ~ s(area_p, bs = 'cr'), data = sub_df_sum, method = 'REML')
p = plot_model(fit, sub_df_sum$area_p, sub_df_sum$W,'area_p')
p_p = p + geom_point(data = sub_df_sum, aes(x = area_p , y = W), color = 'grey14', alpha = 0.3) + 
  labs(x = 'área del papus (mm2)') + scale_x_continuous(breaks = seq(100,250,25)) + theme_bw() 
p_p = p_p + geom_vline(aes(xintercept = mean(sub_df_sum$area_p)), 
                       linetype = 'dashed', color = 'green')
eff2 = wrap_plots(p_m + ggtitle("(a)") +
             p_p + ggtitle("(b)"))


### graficos final #####
pdf('splines_cubicos_indiv.pdf')
eff2
dev.off()
png('splines_cubicos_indiv.png',width = 20, height = 15, units = 'cm', res = 600)
eff2
dev.off()



x_eval <- seq(min(sub_df_sum$masa), max(sub_df_sum$masa), length.out = 1000)

# Predecir valores y derivadas
predicciones <- predict(fit, newdata = data.frame(masa = x_eval), type = "response")

# Usar `deriv` para calcular las derivadas primera y segunda

# Calcular la derivada primera usando diferencias finitas
h <- diff(x_eval)[1] # Asumimos que los puntos de x_eval están equiespaciados
derivada_primera <- diff(predicciones) / h

# Calcular la derivada segunda usando diferencias finitas en la derivada primera
derivada_segunda <- diff(derivada_primera) / h

# Ajustar los vectores para tener la misma longitud que x_eval
derivada_primera <- c(NA, derivada_primera) # Añadir NA al principio para igualar longitud
derivada_segunda <- c(NA, NA, derivada_segunda) # Añadir NAs al principio para igualar longitud

# Encontrar puntos de inflexión donde la derivada segunda cambia de signo
puntos_inflexion <- x_eval[which(diff(sign(derivada_segunda)) != 0) + 2] # Ajustar índice debido a los NAs

# Visualización
plot(sub_df_sum$masa, sub_df_sum$W, main = "Spline cúbica y sus derivadas", xlab = "area_p", ylab = "W")
lines(x_eval, predicciones, col = "blue", lwd = 2, lty = 1)
lines(x_eval, derivada_primera, col = "green", lwd = 2, lty = 2)
lines(x_eval, derivada_segunda, col = "red", lwd = 2, lty = 3)
abline(h = 0, col = "gray", lwd = 0.5)
points(puntos_inflexion, predict(fit, newdata = data.frame(masa = puntos_inflexion)), col = "red", pch = 19)
legend("topright", legend = c("Spline cúbica", "Derivada primera", "Derivada segunda", "Puntos de inflexión"),
       col = c("blue", "green", "red", "red"), lty = c(1, 2, 3, NA), lwd = 2, pch = c(NA, NA, NA, 19))

# Imprimir los puntos de inflexión
puntos_inflexion








x_eval <- seq(min(sub_df_sum$area_p), max(sub_df_sum$area_p), length.out = 1000)

# Predecir valores y derivadas
predicciones <- predict(fit, newdata = data.frame(area_p = x_eval), type = "response")

# Usar `deriv` para calcular las derivadas primera y segunda

# Calcular la derivada primera usando diferencias finitas
h <- diff(x_eval)[1] # Asumimos que los puntos de x_eval están equiespaciados
derivada_primera <- diff(predicciones) / h

# Calcular la derivada segunda usando diferencias finitas en la derivada primera
derivada_segunda <- diff(derivada_primera) / h

# Ajustar los vectores para tener la misma longitud que x_eval
derivada_primera <- c(NA, derivada_primera) # Añadir NA al principio para igualar longitud
derivada_segunda <- c(NA, NA, derivada_segunda) # Añadir NAs al principio para igualar longitud

# Encontrar puntos de inflexión donde la derivada segunda cambia de signo
puntos_inflexion <- x_eval[which(diff(sign(derivada_segunda)) != 0) + 2] # Ajustar índice debido a los NAs

# Visualización
plot(sub_df_sum$area_p, sub_df_sum$W, main = "Spline cúbica y sus derivadas", xlab = "area_p", ylab = "W")
lines(x_eval, predicciones, col = "blue", lwd = 2, lty = 1)
lines(x_eval, derivada_primera, col = "green", lwd = 2, lty = 2)
lines(x_eval, derivada_segunda, col = "red", lwd = 2, lty = 3)
abline(h = 0, col = "gray", lwd = 0.5)
points(puntos_inflexion, predict(fit, newdata = data.frame(area_p = puntos_inflexion)), col = "red", pch = 19)
legend("topright", legend = c("Spline cúbica", "Derivada primera", "Derivada segunda", "Puntos de inflexión"),
       col = c("blue", "green", "red", "red"), lty = c(1, 2, 3, NA), lwd = 2, pch = c(NA, NA, NA, 19))

# Imprimir los puntos de inflexión
puntos_inflexion



### superficie de seleccion gams
## SUPERFICIE DE SELECCION GAMS ####

m3 <- gam(W ~ s(area_p, bs = "cr") + s(masa, bs = "cr" ), 
          data=sub_df_sum, method = "REML", 
          select = TRUE)
m4 <- gam(W ~ te(area_p, masa) , 
          data = sub_df_sum, method = "REML",
          select = TRUE)
m5 <- gam(W ~ s(area_p, bs = "cr") + s(masa, bs = "cr") + ti(area_p, masa),
          data = sub_df_sum, method = "REML", 
          select = TRUE)
m6 = gam(W ~ s(area_p, bs = "cr") + s(masa, bs = "cr") + ti(area_p, masa) ,
         data = sub_df_sum, method = "REML", 
         select = TRUE)
m2 <- gam(W ~ s(area_p, masa, k = 10, bs = "tp") , data=sub_df_sum, 
          method = "REML", 
          select = TRUE)
AIC(m3, m4, m5)
summary(m4)
vis.gam(m2, view = c("masa", "area_p"), type = "response", 
        plot.type = "contour", color = "cm")

new_data = predict_gam(m4 ,length_out = 50)
### graficos final ##### 
ggplot(new_data, aes(x = area_p, y = masa, z = W)) +
  geom_tile(aes(fill = W)) +  # Superficie continua
  geom_contour(aes(z = W),color = 'black', size = 0.3) +  # Contornos blancos
  scale_fill_viridis_c(option = "viridis", name = "Colonización relativa",
                       breaks = seq(0, 2.4, by = 0.2),
                       labels = seq(0, 2.4, by = 0.2)) +  # Escala continua
  labs(x = "Área del Papus (mm2)", 
       y = "Masa de la diáspora (mg)") + scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  theme_bw() + theme(legend.key.height = unit(1.5, "cm"),  # Altura de la barra
                     legend.key.width = unit(0.7, "cm"))

plot_superficie = ggplot(new_data, aes(x = masa, y = area_p, z = W))  + 
  geom_contour_fill(aes(fill = after_stat(level))) + 
  geom_point(data = sub_df_sum, aes(x = masa, y = area_p, z = W), color = 'grey15', alpha = 0.5) +
  scale_fill_discretised(low = "#440154CC", high = "#FDE725CC") +  
  geom_vline(aes(xintercept= mean(sub_df_sum$masa)), 
             linetype = 'dashed', color = 'green') +
  geom_hline(aes(yintercept= mean(sub_df_sum$area_p)), 
             linetype = 'dashed', color = 'green') +
  labs(x = "Masa de la diáspora (mg)", 
       y = "Área del Papus (mm2)",
       fill = "Colonización\n  predicha") + 
  scale_y_continuous(breaks = seq(125, 250, by = 25),
                     limits = c(min(sub_df_sum$area_p), max(sub_df_sum$area_p)),
                     expand = c(0,0)) + 
  scale_x_continuous(breaks = seq(0.5, 2.5, by = .5),
                     limits = c(min(sub_df_sum$masa), max(sub_df_sum$masa)),
                     expand = c(0,0)) +
  theme_bw() + theme(
    legend.key.height = unit(0.8, "cm"),  # Altura de la barra
    legend.key.width = unit(0.5, "cm"),     # Ancho de la barra
    legend.box = "vertical",            # Poner las leyendas en fila
    legend.position = "right"            # Ubicar las leyendas en la parte inferior
  )



pdf('superficie_seleccion_indiv.pdf')
plot_superficie
dev.off()
png('superficie_seleccion_indiv.png',width = 20, height = 15, units = 'cm', res = 600)
plot_superficie
dev.off()
