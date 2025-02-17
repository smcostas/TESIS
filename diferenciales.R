library(tidyverse);library(viridis);library(lme4); library(fitdistrplus); library(MuMIn); library(cAIC4); library(parallel);library(pbkrtest)
library(ggcorrplot); library(PerformanceAnalytics); library(reshape2); library(hier.part); library(patchwork); library(ggeffects)
library(mgcv);library(boot); library(parallel)
## data #####
df <- read.csv('df_sv.csv', header = T)
sub_df <- df[,c('pob','indiv', 'fruto', 'sv', 'prob_germ', 'germ', 'masa', 
                'area_p', 'l_papus', 'vol_c', 'capitulos', 'altura_planta', 'pl_m', 'pl_v')]

df_germ = read.csv('df_gm.csv')

mean(df_germ$masa)
sub_df$inverse_sv <- 1/sub_df$sv

### dispersion ######
sub_df <- sub_df %>%
  group_by(pob, indiv) %>% # Agrupamos por población e individuo
  mutate(indiv2 = cur_group_id()) %>% # Asignamos un ID único continuo
  ungroup()
sub_df_fr <- sub_df %>%
  dplyr::select(pob,indiv,indiv2,fruto, inverse_sv, masa, area_p) %>% 
  group_by(pob,indiv,fruto) %>%
  summarise(across(c( indiv2,inverse_sv, masa, area_p), mean, .names = "{col}"))
sub_df_fr$W = sub_df_fr$inverse_sv/mean(sub_df_fr$inverse_sv)
### germinacion #######
df_germ <- df_germ %>%
  group_by(pob, indiv) %>% # Agrupamos por población e individuo
  mutate(indiv2 = cur_group_id()) %>% # Asignamos un ID único continuo
  ungroup()

## diferenciales lineales para habilidad de dispersion ########

### masa ######

dif_d = lmer(W ~ scale(masa) + (1|indiv2), data = sub_df_fr)
summary(dif_d)
dif_d = lmer(W ~ scale(masa) + (1|pob/indiv), data = sub_df_fr)
summary(dif_d)
#### boostrap para significancia #####

fixef_fun <- function(model) {
  return(fixef(model))
}
get_conf_intervals <- function(results, conf_levels) {
  conf_intervals <- list()
  for (i in 1:ncol(results$t)) {
    ci <- boot.ci(results, index = i, type = "perc", conf = conf_levels)
    conf_intervals[[names(results$t0)[i]]] <- ci ## nombro las sub listas como la variable
  }
  return(conf_intervals)
}
n_cores <- detectCores() - 1
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
boot_results <- bootMer(dif_d, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
conf_intervals <- get_conf_intervals(boot_results, conf_levels) ## re contra siginifcativa

### area p #######

dif_d2 = lmer(W ~ scale(area_p) + (1|indiv2), data = sub_df_fr)
summary(dif_d2)
dif_d2 = lmer(W ~ scale(area_p) + (1|pob/indiv), data = sub_df_fr)
summary(dif_d2)

#### boostrap para significancia #####

boot_results <- bootMer(dif_d2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
conf_intervals <- get_conf_intervals(boot_results, conf_levels) ## re contra siginifcativa

## diferenciales no lineales para habilidad de dispersion ########

### masa ####

dif_d3 = lmer(W ~ I(0.5*(scale(masa)^2)) + (1|indiv2), data = sub_df_fr)
summary(dif_d3)
dif_d3 = lmer(W ~ I(0.5*(scale(masa)^2)) + (1|pob/indiv), data = sub_df_fr)
summary(dif_d3)
#### boot para significancia #####
boot_results <- bootMer(dif_d3, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
conf_intervals <- get_conf_intervals(boot_results, conf_levels) ## re contra siginifcativa

### area p ####
dif_d4 = lmer(W ~ I(0.5*(scale(area_p)^2)) + (1|indiv2), data = sub_df_fr)
summary(dif_d4)
dif_d4 = lmer(W ~ I(0.5*(scale(area_p)^2)) + (1|pob/indiv), data = sub_df_fr)
summary(dif_d4)

#### boot para significancia #####

boot_results <- bootMer(dif_d4, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
conf_intervals <- get_conf_intervals(boot_results, conf_levels) ## no da significativa, tiene sentido


## diferenciales lineales para germinacion ####

dif_g = glmer(germina ~ scale(masa) + (1|indiv2), family = binomial, data = df_germ)
summary(dif_g)

### boostrap para sign #####

boot_results <- bootMer(dif_g, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
conf_intervals <- get_conf_intervals(boot_results, conf_levels) ## no da significativa, tiene sentido


### transformación ######

alpha <- fixef(dif_g)["scale(masa)"]  # Extrae el coeficiente de masa (beta para la masa)
se_alpha <- sqrt(diag(vcov(dif_g)))["scale(masa)"] ## error estándar, podemos comprobar que es igual al que da el summary

W_z <- predict(dif_g, type = "response") ### obtenemos probabilidades estimadas
mean_W <- mean(W_z * (1 - W_z)) ## factor de transformación es la multiplicación de el p(exito)(germina) * 1-p(exito) que es igual al fracaso es decir no germina
beta_avggrad <- mean_W * alpha ## la multiplicación del alpha por la media del W 
se_beta_avggrad <- mean_W * se_alpha ## error estándar para la transformación

## diferenciales no lineales para germinacion ########
dif_g2 = glmer(germina ~ I(0.5*(scale(masa)^2)) + (1|indiv2),family = binomial, data = df_germ)
summary(dif_g2)

### boostrap para sign #####

boot_results2 <- bootMer(dif_g2, fixef_fun, nsim = 1000, use.u = T, re.form = NULL, ncpus = n_cores)
conf_intervals2 <- get_conf_intervals(boot_results2, conf_levels) ## no da significativa, tiene sentido

### transformacion #####
gamma <- fixef(dif_g2)["I(0.5 * (scale(masa)^2))"]
se_gamma <- sqrt(diag(vcov(dif_g2)))["I(0.5 * (scale(masa)^2))"]

W_z <- predict(dif_g2, type = "response")

mean_W2 <- mean((W_z * (1 - W_z))^2) # Calcular el factor de transformación (ahora al cuadrado)
# porque la relación no es lineal en la escala de la regresión logística.

gamma_avggrad <- mean_W2 * gamma
se_gamma_avggrad <- mean_W2 * se_gamma

