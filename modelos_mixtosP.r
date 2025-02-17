library(tidyverse);library(phytools);library(corHMM);library(geiger)
library(viridis); library(nlme); library(ape); library(performance)
library(predictmeans); library(brms); library(DHARMa) ; library(bayesplot)
library(ggeffects); library(patchwork);library(emmeans)

tree <- read.tree('final_tree.tre')
df_sum <- read.csv('final_df.csv', row.names = 1)
summary(df_sum)
which(df_sum$masa==min(df_sum$masa))
which(df_sum$masa==max(df_sum$masa))
df_sum$species = rownames(df_sum)
df_sum[10,]$species
df_sum[46,]$species
which(df_sum$area_p==min(df_sum$area_p))
which(df_sum$area_p==max(df_sum$area_p))
df_sum[34,]$species
df_sum[67,]$species
df_sum[26,]$area_p

which(df_sum$area_p==min(df_sum$area_p))
which(df_sum$area_p==max(df_sum$area_p))
df_sum[34,]$species
df_sum[67,]$species
df_sum[26,]$area_p
sd(df_sum$area_p)
which(df_sum$wind.speed==min(df_sum$wind.speed))
which(df_sum$wind.speed==max(df_sum$wind.speed))
df_sum[54,]$Ecoregion
df_sum[23,]$Ecoregion

df_sum$species = rownames(df_sum)
chk = name.check(tree,df_sum)
df_sum = df_sum[df_sum$species != chk$data_not_tree,]
name.check(tree,df_sum)
df_sum$tipo_papus = factor(df_sum$tipo_papus, levels = c('setose', 'paleaceous', 'aristate', 'epappose'), 
                                                     labels = c('setoso','paleaceo', 'aristado', 'epaposo')) ## correer solo una vez
summary(df_sum$tipo_papus)

### lista de especies #######
df_lista = data.frame(Especies = rownames(df_sum),
                      Tribu = df_sum$Tribu,
                      Tipo_papus = df_sum$tipo_papus)
write.csv(df_lista, 'lista de especies.csv')

##pANOVA ######
A <- ape::vcv.phylo(tree, corr = T)
fit_panova_bayes <- brm(
  sv ~ tipo_papus  + (1 | gr(species, cov = A)),
  data = df_sum,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)

summary(fit_panova_bayes)
r2_bayes(fit_panova_bayes)
emms = emmeans(fit_panova_bayes, ~ tipo_papus)
pares <- pairs(emms)
pares = as.data.frame(pares)
summary_pares <- summary(pares)
letras <- data.frame(tipo_papus = c("aristado", "epaposo", "paleaceo", "setoso"),
                     letras = c("b", "b", "b", "a"))
levels(df_sum$tipo_papus)

colores <- c("aristado" = "#AA337DFF", 
             "epaposo" = "#000004FF", 
             "paleaceo" = "#F7725CFF", 
             "setoso" = "#FDE2A2FF")
n_counts <- df_sum %>%
  group_by(tipo_papus) %>%
  summarise(n = n())

p_boxplot = ggplot(df_sum, aes(x = tipo_papus, y = sv, fill = tipo_papus)) +
  geom_boxplot() +
  scale_fill_manual(values = colores) + 
  geom_text(data = letras, 
            aes(x = tipo_papus, y =c(3.03, 3.03, 2.5, 2.3), label = letras), 
            size = 3, vjust = 0) +  
  geom_text(data = n_counts, 
            aes(x = tipo_papus, y = min(df_sum$sv) - 1e-02, label = paste0("n=", n)), 
            size = 3, vjust = 1.5, color = "black") + 
  labs(x = "Tipo de Papus", y = "Velocidad de asentamiento (m/s)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0),legend.position = "none")

p_boxplot
dif_plot = ggplot(pares, aes(x = contrast, y = estimate, ymin = lower.HPD, ymax = upper.HPD)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "", y = "Diferencia de medias estimada") +
  theme_bw() + coord_flip()
combined_plot <- wrap_plots(p_boxplot, dif_plot, widths = c(1.5,1)) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")")
pdf('panova_brms.pdf',width = 15, height = 8, 'cm')
#png('panova_brms.png',width = 18, height = 10, 'cm', res = 600)
combined_plot
dev.off()

posterior <- as.array(fit_panova_bayes)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_tipo_papuspaleaceo") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "b_tipo_papusaristado") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "b_tipo_papusepaposo") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_tipo_papuspaleaceo", lags = 10)
mcmc_acf(posterior, pars = "b_tipo_papusaristado", lags = 10)
mcmc_acf(posterior, pars = "b_tipo_papusepaposo", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior, pars = c("b_tipo_papuspaleaceo", "b_tipo_papusaristado", 'b_tipo_papusepaposo', 'Intercept'))

pp = posterior_predict(fit_panova_bayes, ndraws = 1000)

y = df_sum$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum$sv)
observed_sd <- sd(df_sum$sv)
observed_max <- max(df_sum$sv)
observed_min <- min(df_sum$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)



qres.fit1=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit1)))

res.m1.brms=cbind(s,
                  df_sum[,c("tipo_papus")],
                  fitted=fitted(fit_panova_bayes, ndraws=1000)[,1],
                  pareto=loo(fit_panova_bayes, pointwise=T)$diagnostics$pareto_k)
colnames(res.m1.brms)[2] = 'tipo_papus'
p1 <- ggplot(res.m1.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m1.brms, aes(x = tipo_papus, y = res)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "tipo de papus", y = "Residuals") +
  theme_minimal()

p4 <- ggplot(res.m1.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()
p5 <- ggplot(res.m1.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2)/(p4+p5)


### extrinsecos * tipo papus diferenciado en setoso vs no setoso #####
df_sum$tipo_papus_b = ifelse(df_sum$tipo_papus == 'setoso', 'setoso', 'no setoso')
df_sum$tipo_papus_b = as.factor(df_sum$tipo_papus_b)
fit_panova_bayes <- brm(
  sv ~ (0 + tipo_papus_b) + (0 + tipo_papus_b:scale(wind.speed))  + (1 | gr(species, cov = A)),
  data = df_sum,
  data2 = list(A = A),  
  family = Gamma(link = "log"),
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)
summary(fit_panova_bayes)

pred_int <- ggpredict(fit_panova_bayes, terms = c("wind.speed [all]", "tipo_papus_b"))
colores_b = c("no setoso" = viridis(1), 
              "setoso" = "#FDE2A2FF")

p_int = ggplot() +
  # Añadir puntos originales del dataset
  # Añadir las bandas de incertidumbre
  geom_ribbon(data = pred_int, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.1) +
  # Añadir las curvas ajustadas
  geom_line(data = pred_int, 
            aes(x = x, y = predicted, color = group), 
            size = 0.9) +
  geom_point(data = df_sum, 
             aes(x = wind.speed, y = sv, color = tipo_papus_b), 
             size = 1.5, alpha = 0.3) +
  scale_color_manual(values = colores_b) +
  scale_fill_manual(values = colores_b) +
  labs(x = "Velocidad media del viento (m/s)", 
       y = "Velocidad de Asentamiento (m/s)", 
       color = "Tipo de papus", 
       fill = "Tipo de papus") +
  theme_minimal() +
  theme(legend.position = "right")
png('modelo interaccion.png',width = 20, height = 15, units = 'cm', res = 600)
p_int
dev.off()
pdf('modelo interaccion.pdf', width = 10, height = 8, 'cm')
p_int
dev.off()

# setosas ############## 
df_sum2 = df_sum %>% filter(tipo_papus == 'setoso')
plot(df_sum2$log.sv ~ df_sum2$air.density)
boxplot(df_sum2$sv)
boxplot(df_sum2$masa)
boxplot(df_sum2$area_p)
boxplot(sqrt(df_sum2$pl_m))
df_sum2 = df_sum2[which(df_sum2$masa < 5 & df_sum2$pl_m < 0.09  & df_sum2$area_p < 750),]

which(df_sum2$wind.speed==min(df_sum2$wind.speed))
which(df_sum2$wind.speed==max(df_sum2$wind.speed))
df_sum2[38,]$wind.speed
df_sum2[12,]$wind.speed
sd(df_sum[,]$wind.speed)
which(df_sum$bio12==min(df_sum$bio12))
which(df_sum$bio12==max(df_sum$bio12))
df_sum[c(40,48),]$bio12
df_sum[9,]$bio12
sd(df_sum[,]$bio12)
which(df_sum$cover==min(df_sum$cover))
which(df_sum$cover==max(df_sum$cover))
df_sum[which(df_sum$cover==min(df_sum$cover)),]$Ecoregion
df_sum[54,]$Ecoregion
sd(df_sum[,]$Ecoregion)
df_sum[which(df_sum$cover==min(df_sum$cover)),]$cod_gps
which(df_sum$air.density==min(df_sum$air.density))
which(df_sum$air.density==max(df_sum$air.density))
df_sum[which(df_sum$air.density==max(df_sum$air.density)),]$Ecoregion
df_sum[which(df_sum$air.density==min(df_sum$air.density)),]$Ecoregion

chk = name.check(tree,df_sum2)
tree2 = drop.tip(tree,chk$tree_not_data)
#df_sum2 = df_sum2[df_sum2$species != chk$data_not_tree,]
name.check(tree2,df_sum2)
summary(df_sum2)

boxplot(df_sum2$sv)
boxplot(df_sum2$masa)
boxplot(df_sum2$area_p)
boxplot(df_sum2$pl_m)
cor(df_sum2$masa,df_sum2$area_p)
cor(df_sum2$masa, sqrt(df_sum2$pl_m))
cor(df_sum2$area_p, df_sum2$pl_m)

#binary_tree <- multi2di(tree2)## uso el binario
A <- ape::vcv.phylo(tree2, corr = T) 

## MODELO 1 #######################################################################

get_prior(form = sv ~ scale(masa) +scale(area_p) + (1 | gr(species, cov = A)),
          data = df_sum2,data2 = list(A = A), family = Gamma(link = "log")) 

prior.m1 = c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("normal(log(0.6), 2)", class = "Intercept")
)

fit1_bayes <- brm(
  sv ~ scale(masa) + scale(area_p)  + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  prior = prior.m1,
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)

summary(fit1_bayes)
r2_bayes(fit1_bayes)
preds_area=ggpredict(fit1_bayes,type="fixed", terms =
                       c("area_p[all]"))
p_ap = ggplot(preds_area, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "area del papus (mm2)", y = "") + geom_point(data = df_sum2, aes(x = area_p, y = sv)) + theme_bw()
p_ap

preds = ggpredict(fit1_bayes, c('masa[all]'))
p_m = ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "masa (mg)", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = masa, y = sv)) + theme_bw()

p_ap + p_m
color_scheme_set('purple')
pi_ma = mcmc_plot(fit1_bayes, prob=0.5, prob_outer = 0.95, regex_pars = c("^b")) + 
  geom_vline(aes(xintercept  = 0.0), linetype = 'solid', color = viridis(1, begin = 0.5)) + 
  scale_y_discrete(labels = c('Intercepto', 'b_ masa', 'b_ área papus')) + 
  theme_bw() + 
  theme(axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())
pi_ma 

posterior <- as.array(fit1_bayes)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalearea_p") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "b_scalemasa") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalemasa", lags = 10)
mcmc_acf(posterior, pars = "b_scalearea_p", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

posterior[1,1,]
dim(posterior)

mcmc_dens_overlay(posterior, pars = c("b_scalemasa", "b_scalearea_p", 'Intercept'))

pp = posterior_predict(fit1_bayes, ndraws = 1000)

y = df_sum2$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum2$sv)
observed_sd <- sd(df_sum2$sv)
observed_max <- max(df_sum2$sv)
observed_min <- min(df_sum2$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)



qres.fit1=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum2$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit1)))

res.m1.brms=cbind(s,
                  df_sum2[,c("masa","area_p")],
                  fitted=fitted(fit1_bayes, ndraws=1000)[,1],
                  pareto=loo(fit1_bayes, pointwise=T)$diagnostics$pareto_k)
p1 <- ggplot(res.m1.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m1.brms, aes(x = scale(masa), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Masa (std)", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m1.brms, aes(x = scale(area_p), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Area_p (std)", y = "Residuals") +
  theme_minimal()
p4 <- ggplot(res.m1.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()
p5 <- ggplot(res.m1.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2)/(p3+p4+p5)
 
## MODELO 2 ... without pareto k extreme points #######################################
loo_fit1 <- loo(fit1_bayes)
pareto_k <- loo_fit1$diagnostics$pareto_k
problematic_points <- which(pareto_k > 0.70)
spp = df_sum2[problematic_points,]$species
df_sum_f_pm = df_sum2[!df_sum2$species%in%spp, ]
#df_sum_f_pm = df_sum_f_pm[which(df_sum_f_pm$sv < 1),]

chk = name.check(tree2,df_sum_f_pm)
tree_f_pm = drop.tip(tree2,chk$tree_not_data)
name.check(tree_f_pm,df_sum_f_pm)

#binary_tree <- multi2di(tree2)
A <- ape::vcv.phylo(tree_f_pm, corr = T)



fit1_bayes_loo <- brm(
  sv ~ scale(masa) + scale(area_p) + (1 | gr(species, cov = A)),
  data = df_sum_f_pm,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)
summary(fit1_bayes_loo)
r2_bayes(fit1_bayes_loo)
preds_area_loo=ggpredict(fit1_bayes_loo,type="fixed", terms =
                       c("area_p[all]"))
p_ap = ggplot(preds_area_loo, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "area del papus (mm2)", y = "") + geom_point(data = df_sum_f_pm, aes(x = area_p, y = sv)) + theme_bw()
p_ap

preds_loo = ggpredict(fit1_bayes_loo, c('masa[all]'))
p_m = ggplot(preds_loo, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "masa (mg)", y = "Velocidad de Asentamiento (m/s)") + geom_point(data = df_sum_f_pm, aes(x = masa, y = sv)) + theme_bw()

p_ap + p_m

posterior_loo <- as.array(fit1_bayes_loo)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior_loo, pars = "b_scalearea_p") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior_loo, pars = "b_scalemasa") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior_loo, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior_loo, pars = "b_scalemasa", lags = 10)
mcmc_acf(posterior_loo, pars = "b_scalearea_p", lags = 10)
mcmc_acf(posterior_loo, pars = "Intercept", lags = 10)
mcmc_acf(posterior_loo, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior_loo, pars = c("b_scalemasa", "b_scalearea_p", 'Intercept'))

pp = posterior_predict(fit1_bayes_loo, ndraws = 1000)

y = df_sum_f_pm$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum_f_pm$sv)
observed_sd <- sd(df_sum_f_pm$sv)
observed_max <- max(df_sum_f_pm$sv)
observed_min <- min(df_sum_f_pm$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)



qres.fit2=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum_f_pm$sv,
                       fittedPredictedResponse=apply(pp, 2, median))
length(qres.fit2)
s2= data.frame(res = qnorm(residuals(qres.fit2)))

res.m2.brms=cbind(s2,
                  df_sum_f_pm[,c("masa","area_p")],
                  fitted=fitted(fit1_bayes_loo, ndraws=1000)[,1],
                  pareto=loo(fit1_bayes_loo, pointwise=T)$diagnostics$pareto_k)

p1 <- ggplot(res.m2.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m2.brms, aes(x = scale(masa), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Masa (std)", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m2.brms, aes(x = scale(area_p), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Area_p (std)", y = "Residuals") +
  theme_minimal()
p4 <- ggplot(res.m2.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()
p5 <- ggplot(res.m2.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2)/(p3+p4+p5)

### graficos final  ######
### area_p ####
preds_area_loo=ggpredict(fit1_bayes_loo,type="fixed", terms =
                           c("area_p[25,50,100,150,200,250,300,350,400,450,500,550,575]"))
preds_area_loo$model = 'Modelo sin puntos influyentes'
preds_area$model = 'Modelo completo'
p_ap = ggplot(preds_area_loo, aes(x = x, y = predicted)) +
  geom_line(color = 'grey35', linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = 'grey35', alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Área del papus (mm2)", y = "Velocidad de asentamiento (m/s)") 
p_ap = p_ap  + geom_line(data = preds_area, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds_area, aes(x = x, y = predicted,ymin = conf.low, ymax = conf.high),fill = viridis(1), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = area_p, y = sv), alpha = 0.2) +
  theme_bw() + theme(axis.title.y =  element_blank())
p_ap
### masa ####
preds_loo=ggpredict(fit1_bayes_loo,type="fixed", terms =
                       c("masa[0.05,0.1,0.2,0.4,0.6,.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.1,3.2,3.3,3.4,3.5,3.7,3.9]"))
p = ggplot(preds_loo, aes(x = x, y = predicted)) +
  geom_line(color = "grey35", linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = "grey35") +  # Banda de intervalo de confianza
  labs(x = "Masa de la diáspora (mg)", y = "Velocidad de asentamiento (m/s)") 
p_m = p  + geom_line(data = preds, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = viridis(1)) + 
  geom_point(data = df_sum2, aes(x = masa, y = sv), alpha = 0.2) + theme_bw()

p_m + p_ap
pdf('model_masa_area.pdf',width = 15, height = 10, 'cm')
(p_m + p_ap)/p_pl
dev.off()

## MODELO pl ######
A <- ape::vcv.phylo(tree2, corr = T)
df_sum2$sqrt_pl = sqrt(df_sum2$pl_m)

get_prior(form = sv ~ scale(sqrt_pl) + (1 | gr(species, cov = A)),
          data = df_sum2,data2 = list(A = A), family = Gamma(link = "log")) 

prior.m2 = c(set_prior("normal(0,2)", class="b", coef = 'scalesqrt_pl'),
             set_prior("normal(log(0.6),2)", class = "Intercept"),
             set_prior("gamma(0.01,0.01)", class = "shape"))

#c(
#prior(student_t(3, 0, 10), class = "b"),  # Priors on regression coefficients
#prior(inv_gamma(2, 0.1), class = "sd")    # Inverse gamma prior on standard deviation
#)

plot(sv~sqrt_pl, data = df_sum2)
boxplot(df_sum2$sqrt_pl)
fit3_bayes <- brm(
  sv ~ scale(sqrt_pl) + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),  # Gamma family for positive continuous data
  prior = prior.m2 ,
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(adapt_delta = 0.999),
  cores = 11
)
summary(fit3_bayes)
r2_bayes(fit3_bayes)
pp_check(fit3_bayes)
posterior <- as.array(fit3_bayes)
posterior[1,,]

preds_pl=ggpredict(fit3_bayes,type="fixed", terms =
                     c("sqrt_pl[all]"))
p_pl = ggplot(preds_pl, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "raiz cuadrada de PL", y = "Velocidad de asentamiento (m/s)") + geom_point(data = df_sum2, aes(x = sqrt_pl, y = sv)) + theme_bw()
p_pl

color_scheme_set('purple')
pi_pl = mcmc_plot(fit3_bayes, prob=0.5, prob_outer = 0.95, regex_pars = c("^b")) + 
  geom_vline(aes(xintercept  = 0.0), linetype = 'solid', color = viridis(1, begin = 0.5)) +  
  scale_y_discrete(labels = c('Intercepto', 'b_ PL')) + 
  theme_bw() + 
  theme(axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())
pi_pl
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalesqrt_pl") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalesqrt_pl", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

posterior[1,1,]
dim(posterior)

mcmc_dens_overlay(posterior, pars = c("b_scalesqrt_pl", 'Intercept'))


pp = posterior_predict(fit3_bayes, ndraws = 1000)
dim(pp)
y = df_sum2$sv
ppc_dens_overlay(y, pp)

pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum2$sv)
observed_sd <- sd(df_sum2$sv)
observed_max <- max(df_sum2$sv)
observed_min <- min(df_sum2$sv)




# Graficar la media de las simulaciones y marcar la media observada
p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)
## el modelo tiene grandes problemas con el maximo

qres.fit3=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum2$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit3)))

res.m3.brms=cbind(s,
                  df_sum2$sqrt_pl,
                  fitted=fitted(fit3_bayes, ndraws=1000)[,1],
                  pareto=loo(fit3_bayes, pointwise=T)$diagnostics$pareto_k)
colnames(res.m3.brms)[2] = 'sqrt_pl_m'
p1 <- ggplot(res.m3.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m3.brms, aes(x = scale(sqrt_pl_m), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "pl", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m3.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()

p4 <- ggplot(res.m3.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2)/(p3+p4)

### paraeto k #####
loo_fit3 <- loo(fit3_bayes)
loo_compare(loo_fit1, loo_fit3) 

pareto_k <- loo_fit3$diagnostics$pareto_k
problematic_points <- which(pareto_k > 0.70)
spp = df_sum2[problematic_points,]$species
df_sum_f_pl = df_sum2[!df_sum2$species%in%spp, ]


## MODELO 4 pl pareto k #####

chk = name.check(tree2,df_sum_f_pl)
tree_f_pl = drop.tip(tree2,chk$tree_not_data)
name.check(tree_f_pl,df_sum_f_pl)

#binary_tree <- multi2di(tree2)
A <- ape::vcv.phylo(tree_f_pl, corr = T)

fit3_bayes_loo <- brm(
  sv ~ scale(sqrt_pl) + (1 | gr(species, cov = A)),
  data = df_sum_f_pl,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),  # Gamma family for positive continuous data
  prior = prior.m2 ,
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(adapt_delta = 0.999),
  cores = 11
)
summary(fit3_bayes_loo)
r2_bayes(fit3_bayes_loo)
posterior <- as.array(fit3_bayes_loo)

preds_pl_loo=ggpredict(fit3_bayes_loo,type="fixed", terms =
                     c("sqrt_pl[all]"))
p_pl = ggplot(preds_pl_loo, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "raiz cuadrada de PL", y = "Velocidad de asentamiento (m/s)") + geom_point(data = df_sum_f_pl, aes(x = sqrt_pl, y = sv)) + theme_bw()
p_pl
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalesqrt_pl") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalesqrt_pl", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior, pars = c("b_scalesqrt_pl", 'Intercept'))


pp = posterior_predict(fit3_bayes_loo, ndraws = 1000)
y = df_sum_f_pl$sv
ppc_dens_overlay(y, pp)

pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum_f_pl$sv)
observed_sd <- sd(df_sum_f_pl$sv)
observed_max <- max(df_sum_f_pl$sv)
observed_min <- min(df_sum_f_pl$sv)




# Graficar la media de las simulaciones y marcar la media observada
p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)
## el modelo tiene grandes problemas con el maximo

qres.fit4=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum_f_pl$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit4)))

res.m4.brms=cbind(s,
                  df_sum_f_pl$sqrt_pl,
                  fitted=fitted(fit3_bayes_loo, ndraws=1000)[,1],
                  pareto=loo(fit3_bayes_loo, pointwise=T)$diagnostics$pareto_k)
colnames(res.m4.brms)[2] = 'sqrt_pl_m'
p1 <- ggplot(res.m4.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m4.brms, aes(x = scale(sqrt_pl_m), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "pl", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m4.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()

p4 <- ggplot(res.m4.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2)/(p3+p4)

### graficos final  ######
### sqrt_PL 
min(df_sum2$sqrt_pl)
preds_pl_loo=ggpredict(fit3_bayes_loo,type="fixed", terms =
                         c("sqrt_pl[all]"))
preds_pl_loo$model = 'modelo pareto K'
p_pl = ggplot(preds_pl_loo, aes(x = x, y = predicted)) +
  geom_line(aes(color = model, linetype = model)) +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Raíz cuadrado de PL", y = "Velocidad de asentamiento (m/s)") 
p_pl
preds_pl=ggpredict(fit3_bayes,type="fixed", terms =
                         c("sqrt_pl[all]"))
preds_pl$model = 'modelo completo'
p_pl = p_pl + geom_line(data = preds_pl, aes(x = x, y = predicted, color = model,linetype = model)) + 
  geom_ribbon(data = preds_pl, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = sqrt_pl, y = sv), alpha = 0.2) +
  scale_color_manual(values = c(viridis(1),'grey35')) + 
  scale_fill_manual(values = c(viridis(1),'grey35')) + 
  scale_linetype_manual(values = c('solid','dashed'))+
  theme_bw() + 
  theme(axis.title.y =  element_blank(),
        legend.title = element_blank(), legend.position = c(0,1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.justification = c(0,1))

combined_plot <- (p_m + p_ap + p_pl) / (pi_ma + pi_pl ) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")") 

combined_plot <- wrap_plots(wrap_plots(
  p_m, p_ap, p_pl,ncol =3),               # Primera fila con 3 gráficos individuales
  wrap_plots(pi_ma, pi_pl, ncol=2),  # Segunda fila combinada en 1 gráfico que ocupa toda la fila
  heights = c(2.5, 1)    # Distribución de alturas: primera fila más alta
) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")")

combined_plot

pdf('models_intrinsecos.pdf',width = 10, height = 8, 'cm')
combined_plot
dev.off()
png('models_intrinsecos.png',width = 20, height = 15, units = 'cm', res = 600)
combined_plot
dev.off()
combined_plot


# CARACTERES EXTRINSECOS ##########################################################
## MODELO 1 ######################################################################
A <- ape::vcv.phylo(tree2, corr = T)
get_prior(form = sv ~  scale(wind.speed) + scale(bio12) + scale(cover) + (1 | gr(species, cov = A)),
          data = df_sum2, data2 = list(A = A), family = Gamma(link = "log"))
prior.m7 = c(set_prior("cauchy(0,2)", class="b"),
             #set_prior("normal(0,2)", class="b", coef = 'scalemasa:scalearea_p'),
             set_prior("cauchy(log(0.5),2)", class = "Intercept"),
             set_prior("gamma(0.01,0.01)", class = "shape"))

fit7_bayes <- brm(
  sv ~  scale(wind.speed) + scale(bio12) + scale(cover) + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),  # Gamma family for positive continuous data
  prior = prior.m7 ,
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(adapt_delta = 0.99),
  cores = 11
)

fit8_bayes <- brm(
  sv ~  scale(wind.speed) + scale(bio12) + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),  # Gamma family for positive continuous data
  prior = prior.m7 ,
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(adapt_delta = 0.99),
  cores = 11
)
loo_fit7 <- loo(fit7_bayes)
loo_fit8 <- loo(fit8_bayes)
loo_compare(loo_fit7, loo_fit8) 
fit7_bayes <- brm(
  sv ~  scale(wind.speed) + scale(bio12) + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),  # Gamma family for positive continuous data
  prior = prior.m7 ,
  warmup=3000,chains=3, iter=10000, thin=5,
  control = list(adapt_delta = 0.99),
  cores = 11
)

summary(fit7_bayes)
bayes_R2(fit7_bayes, summary = T)
r2_bayes(fit7_bayes)
library(ggeffects)
preds_w=ggpredict(fit7_bayes,type="fixed", terms =
                    c("wind.speed[all]"))
p_w = ggplot(preds_w, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "velocidad del viento (m/s)", y = "VA (m/s)") + geom_point(data = df_sum2, aes(x = wind.speed, y = sv)) + theme_bw()
p_w
preds_p=ggpredict(fit7_bayes,type="fixed", terms =
                    c("bio12[all]"))
p_p = ggplot(preds_p, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Presipitación media anual (mm)", y = "VA (m/s)") + geom_point(data = df_sum2, aes(x = bio12, y = sv)) + theme_bw()
p_w + p_p
color_scheme_set('purple')
pi_ext = mcmc_plot(fit7_bayes, prob=0.5, prob_outer = 0.95, regex_pars = c("^b")) + 
  geom_vline(aes(xintercept  = 0.0), linetype = 'solid', color = viridis(1, begin = 0.5)) +  
  scale_y_discrete(labels = c('Intercepto', 'b_Vel viento', 'b_Precipitación MA')) + 
  theme_bw() + 
  theme(axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())
pi_ext

posterior <- as.array(fit7_bayes)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalewind.speed") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalesqrt_pl", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior, pars = c("b_scalebio12", 'b_scalewind.speed','Intercept'))

pp = posterior_predict(fit7_bayes, ndraws = 500)

y = df_sum2$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum2$sv)
observed_sd <- sd(df_sum2$sv)
observed_max <- max(df_sum2$sv)
observed_min <- min(df_sum2$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)





qres.fit7=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum2$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit7)))

res.m7.brms=cbind(s,
                  df_sum2[,c("wind.speed","bio12")],
                  fitted=fitted(fit7_bayes, ndraws=1000)[,1],
                  pareto=loo(fit7_bayes, pointwise=T)$diagnostics$pareto_k)
p1 <- ggplot(res.m7.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m7.brms, aes(x = scale(wind.speed), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "velocidad de viento", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m7.brms, aes(x = scale(bio12), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "presipitación media anual", y = "Residuals") +
  theme_minimal()

#p4 <- ggplot(res.m7.brms, aes(x = scale(cover), y = res)) +
#  geom_point() +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#  labs(x = "cobertura de árboles", y = "Residuals") +
#  theme_minimal()

p5 <- ggplot(res.m7.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()

p6 <- ggplot(res.m7.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2+p3)/(p5+p6)

## MODEL 2 pareto k ######
loo_fit7 <- loo(fit7_bayes)
pareto_k <- loo_fit7$diagnostics$pareto_k
problematic_points <- which(pareto_k > 0.70)
spp = df_sum2[problematic_points,]$species
df_sum_f_ev = df_sum2[!df_sum2$species%in%spp, ]


chk = name.check(tree2,df_sum_f_ev)
tree_f_ev = drop.tip(tree2,chk$tree_not_data)
name.check(tree_f_ev,df_sum_f_ev)

#binary_tree <- multi2di(tree2)
A <- ape::vcv.phylo(tree_f_ev, corr = T)

fit7_bayes_loo <- brm(
  sv ~  scale(wind.speed) + scale(bio12) + (1 | gr(species, cov = A)),
  data = df_sum_f_ev,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),  # Gamma family for positive continuous data
  prior = prior.m7 ,
  warmup=5000,chains=3, iter=15000, thin=5,
  control = list(adapt_delta = 0.99),
  cores = 11
)

summary(fit7_bayes_loo)
bayes_R2(fit7_bayes_loo, summary = T)
r2_bayes(fit7_bayes_loo)
library(ggeffects)
preds_w_loo=ggpredict(fit7_bayes_loo,type="fixed", terms =
                    c("wind.speed[all]"))
p_w = ggplot(preds_w_loo, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "velocidad del viento (m/s)", y = "VA (m/s)") + geom_point(data = df_sum_f_ev, aes(x = wind.speed, y = sv)) + theme_bw()
p_w
preds_p=ggpredict(fit7_bayes_loo,type="fixed", terms =
                    c("bio12[all]"))
p_p = ggplot(preds_p, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Presipitación media anual", y = "VA (m/s)") + geom_point(data = df_sum_f_ev, aes(x = bio12, y = sv)) + theme_bw()
p_w + p_p


posterior <- as.array(fit7_bayes_loo)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalewind.speed") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalesqrt_pl", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior, pars = c("b_scalebio12", 'b_scalewind.speed','Intercept'))

pp = posterior_predict(fit7_bayes_loo, ndraws = 500)

y = df_sum_f_ev$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum_f_ev$sv)
observed_sd <- sd(df_sum_f_ev$sv)
observed_max <- max(df_sum_f_ev$sv)
observed_min <- min(df_sum_f_ev$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", size = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)





qres.fit7=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum_f_ev$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit7)))

res.m7.brms=cbind(s,
                  df_sum_f_ev[,c("wind.speed","bio12")],
                  fitted=fitted(fit7_bayes_loo, ndraws=1000)[,1],
                  pareto=loo(fit7_bayes_loo, pointwise=T)$diagnostics$pareto_k)
p1 <- ggplot(res.m7.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m7.brms, aes(x = scale(wind.speed), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "velocidad de viento", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m7.brms, aes(x = scale(bio12), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "presipitación media anual", y = "Residuals") +
  theme_minimal()

#p4 <- ggplot(res.m7.brms, aes(x = scale(cover), y = res)) +
#  geom_point() +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#  labs(x = "cobertura de árboles", y = "Residuals") +
#  theme_minimal()

p5 <- ggplot(res.m7.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()

p6 <- ggplot(res.m7.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2+p3)/(p5+p6)

## GRAFICOS FINALES EXTRINSECOS #########################################################################################
preds_ws_loo=ggpredict(fit7_bayes_loo,type="fixed", terms =
                         c("wind.speed[all]"))
preds_ws_loo$model = 'modelo Pareto K'
p_ws = ggplot(preds_ws_loo, aes(x = x, y = predicted)) +
  geom_line(aes(color = model, linetype = model)) +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Velocidad del viento (m/s)", y = "Velocidad de asentamiento (m/s)") 
p_ws
preds_ws=ggpredict(fit7_bayes,type="fixed", terms =
                     c("wind.speed[all]"))
preds_ws$model = 'modelo completo'
p_ws = p_ws + geom_line(data = preds_ws, aes(x = x, y = predicted, color = model,linetype = model)) + 
  geom_ribbon(data = preds_ws, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = wind.speed, y = sv), alpha = 0.2) +
  scale_color_manual(values = c(viridis(1),'grey35')) + 
  scale_fill_manual(values = c(viridis(1),'grey35')) + 
  scale_linetype_manual(values = c('solid','dashed'))+
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(), legend.position = c(0,1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.justification = c(0,1))

preds_pma_loo=ggpredict(fit7_bayes_loo,type="fixed", terms =
                         c("bio12[all]"))
preds_pma_loo$model = 'modelo Pareto K'
p_pma = ggplot(preds_pma_loo, aes(x = x, y = predicted)) +
  geom_line(color = 'grey35',linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = 'grey35', alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Precipitación media anual (mm)", y = "Velocidad de asentamiento (m/s)") 

preds_pma=ggpredict(fit7_bayes,type="fixed", terms =
                     c("bio12[all]"))
preds_pma$model = 'modelo completo'
p_pma = p_pma + geom_line(data = preds_pma, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds_pma, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),fill = viridis(1), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = bio12, y = sv), alpha = 0.2) +
  theme_bw()

combined_plot <- wrap_plots(wrap_plots(p_pma,p_ws,ncol = 2), wrap_plots(plot_spacer(),pi_ext,plot_spacer(), ncol = 3, widths = c(1,3,1) ), nrow = 2, heights = c(1.5,1)) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")")

#pdf('models_extrinsecos.pdf',width = 15, height = 10, 'cm')

png('models_extrinsecos.png',width = 22, height = 18, 'cm', res =600)

combined_plot
dev.off()

## MODELO ALL TOGETHER ########


get_prior(form =sv ~  scale(masa) + scale(area_p) + scale(wind.speed) + scale(bio12)  + (1 | gr(species, cov = A)),
          data = df_sum2,data2 = list(A = A), family = Gamma(link = "log")) 

prior.m9 = c(
  set_prior("cauchy(0, 0.5)", class = "b"),
  set_prior("normal(log(0.6), 2)", class = "Intercept")
)

fit8_bayes <- brm(
  sv ~  scale(masa) + scale(area_p) + scale(wind.speed) + scale(bio12)  + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  prior = prior.m9,
  warmup=10000,chains=3, iter=30000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)
summary(fit8_bayes)
r2_bayes(fit8_bayes)

preds_w=ggpredict(fit8_bayes,type="fixed", terms =
                    c("wind.speed[all]"))
p_w = ggplot(preds_w, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "velocidad del viento (m/s)", y = "VA (m/s)") + geom_point(data = df_sum2, aes(x = wind.speed, y = sv)) + theme_bw()
p_w
preds_p=ggpredict(fit8_bayes,type="fixed", terms =
                    c("bio12[all]"))
p_p = ggplot(preds_p, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +  # Banda de intervalo de confianza
  labs(x = "Presipitación media anual (mm)", y = "VA (m/s)") + geom_point(data = df_sum2, aes(x = bio12, y = sv)) + theme_bw()
p_w + p_p





posterior <- as.array(fit8_bayes)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalearea_p") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "b_scalemasa") + 
  xlab("Post-warmup iteration") 
mcmc_trace(posterior, pars = "b_scalewind.speed") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "b_scalebio12") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalemasa", lags = 10)
mcmc_acf(posterior, pars = "b_scalearea_p", lags = 10)
mcmc_acf(posterior, pars = "b_scalewind.speed", lags = 10)
mcmc_acf(posterior, pars = "b_scalebio12", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior, pars = c("b_scalemasa", "b_scalearea_p", 'b_scalewind.speed',
                                       'b_scalebio12', 'Intercept','shape'))

pp = posterior_predict(fit8_bayes, ndraws = 1000)

y = df_sum2$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum2$sv)
observed_sd <- sd(df_sum2$sv)
observed_max <- max(df_sum2$sv)
observed_min <- min(df_sum2$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)



qres.fit8=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum2$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit8)))

res.m8.brms=cbind(s,
                  df_sum2[,c("masa","area_p", "wind.speed", "bio12")],
                  fitted=fitted(fit8_bayes, ndraws=1000)[,1],
                  pareto=loo(fit8_bayes, pointwise=T)$diagnostics$pareto_k)
p1 <- ggplot(res.m8.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m8.brms, aes(x = scale(masa), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Masa (std)", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m8.brms, aes(x = scale(area_p), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Area_p (std)", y = "Residuals") +
  theme_minimal()
p4 <- ggplot(res.m8.brms, aes(x = scale(wind.speed), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "wind speed (std)", y = "Residuals") +
  theme_minimal()
p5 <- ggplot(res.m8.brms, aes(x = scale(bio12), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "temperatura (std)", y = "Residuals") +
  theme_minimal()
p6 <- ggplot(res.m8.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()
p7 <- ggplot(res.m8.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2+p3)/(p4+p5+p6+p7)

## MODEL 8 PARETOK #######
loo_fit8 <- loo(fit8_bayes)
pareto_k <- loo_fit8$diagnostics$pareto_k
problematic_points <- which(pareto_k > 0.70)
spp = df_sum2[problematic_points,]$species
df_sum_f_at = df_sum2[!df_sum2$species%in%spp, ]
chk = name.check(tree2,df_sum_f_at)
tree_at = drop.tip(tree2,chk$tree_not_data)
name.check(tree_at,df_sum_f_at)

#binary_tree <- multi2di(tree2)
A <- ape::vcv.phylo(tree_at, corr = T)
fit8_bayes_loo <- brm(
  sv ~  scale(masa) + scale(area_p) + scale(wind.speed) + scale(bio12)  + (1 | gr(species, cov = A)),
  data = df_sum_f_at,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  prior = prior.m9,
  warmup=10000,chains=3, iter=20000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)
summary(fit8_bayes_loo)


preds_area= ggpredict(fit8_bayes,type="fixed", terms =
            c("area_p[25,50,100,150,200,250,300,350,400,450,500,550,575]"))

preds_area_loo=ggpredict(fit8_bayes_loo,type="fixed", terms =
                           c("area_p[25,50,100,150,200,250,300,350,400,450,500,550,575]"))
preds_area_loo$model = 'Modelo sin puntos influyentes'
preds_area$model = 'Modelo completo'
p_ap = ggplot(preds_area_loo, aes(x = x, y = predicted)) +
  geom_line(color = 'grey35', linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = 'grey35', alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Área del papus (mm2)", y = "Velocidad de asentamiento (m/s)") 
p_ap = p_ap  + geom_line(data = preds_area, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds_area, aes(x = x, y = predicted,ymin = conf.low, ymax = conf.high),fill = viridis(1), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = area_p, y = sv), alpha = 0.2) +
  theme_bw() + theme(axis.title.y =  element_blank())
p_ap
### graficos final ####
preds = ggpredict(fit8_bayes,type="fixed", terms =
                    c("masa[0.05,0.1,0.2,0.4,0.6,.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.1,3.2,3.3,3.4,3.5,3.7,3.9]"))
preds_loo=ggpredict(fit8_bayes_loo,type="fixed", terms =
                      c("masa[0.05,0.1,0.2,0.4,0.6,.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.1,3.2,3.3,3.4,3.5,3.7,3.9]"))
p = ggplot(preds_loo, aes(x = x, y = predicted)) +
  geom_line(color = "grey35", linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = "grey35") +  # Banda de intervalo de confianza
  labs(x = "Masa de la diáspora (mg)", y = "Velocidad de asentamiento (m/s)") 
p_m = p  + geom_line(data = preds, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = viridis(1)) + 
  geom_point(data = df_sum2, aes(x = masa, y = sv), alpha = 0.2) + theme_bw()

preds_wind = ggpredict(fit8_bayes,type="fixed", terms =
                    c("wind.speed[all]"))
preds_wind_loo=ggpredict(fit8_bayes_loo,type="fixed", terms =
                      c("wind.speed[all]"))
p = ggplot(preds_wind_loo, aes(x = x, y = predicted)) +
  geom_line(color = "grey35", linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = "grey35") +  # Banda de intervalo de confianza
  labs(x = "Velocidad media del viento (m/s)", y = "Velocidad de asentamiento (m/s)") 
p_w = p  + geom_line(data = preds_wind, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds_wind, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = viridis(1)) + 
  geom_point(data = df_sum2, aes(x = wind.speed, y = sv), alpha = 0.2) + theme_bw() +  theme(axis.title.y =  element_blank())


preds_t_loo=ggpredict(fit8_bayes_loo,type="fixed", terms =
                         c("bio12[all]"))
preds_t_loo$model = 'modelo pareto K'
p_t = ggplot(preds_t_loo, aes(x = x, y = predicted)) +
  geom_line(aes(color = model, linetype = model)) +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Precipitación media anual (mm)", y = "Velocidad de asentamiento (m/s)") 
p_t
preds_t=ggpredict(fit8_bayes,type="fixed", terms =
                     c("bio12[all]"))
preds_t$model = 'modelo completo'
p_t = p_t + geom_line(data = preds_t, aes(x = x, y = predicted, color = model,linetype = model)) + 
  geom_ribbon(data = preds_t, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = bio12, y = sv), alpha = 0.2) +
  scale_color_manual(values = c(viridis(1),'grey35')) + 
  scale_fill_manual(values = c(viridis(1),'grey35')) + 
  scale_linetype_manual(values = c('solid','dashed'))+
  theme_bw() + 
  theme(axis.title.y =  element_blank(),
        legend.title = element_blank(), legend.position = c(1,1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.justification = c(1,1))
p_t
color_scheme_set('purple')
pi = mcmc_plot(fit8_bayes, prob=0.5, prob_outer = 0.95, regex_pars = c("^b")) + 
  geom_vline(aes(xintercept  = 0.0), linetype = 'solid', color = viridis(1, begin = 0.5)) + 
  scale_y_discrete(labels = c('intercepto', 'masa', 'área papus', 'vel viento', 'precipitación MA')) + 
  theme_bw() + 
  theme(axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())
pi 




combined_plot <- wrap_plots(wrap_plots(
  p_m, p_ap, p_w, p_t,ncol =4),               # Primera fila con 3 gráficos individuales
  wrap_plots(plot_spacer(),pi, plot_spacer(), ncol=3),  # Segunda fila combinada en 1 gráfico que ocupa toda la fila
  heights = c(2.5, 1)    # Distribución de alturas: primera fila más alta
) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")")

pdf('modelo completo.pdf',width = 10, height = 8, 'cm')
combined_plot
dev.off()

png('modelo completo.png',width = 30, height = 15, units = 'cm', res = 600)

combined_plot
dev.off()


## sv ~ pl + extrinsecas ######
fit9_bayes <- brm(
  sv ~  scale(sqrt_pl) + scale(wind.speed) + scale(bio12)  + (1 | gr(species, cov = A)),
  data = df_sum2,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  prior = prior.m9,
  warmup=10000,chains=3, iter=20000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)
summary(fit9_bayes)
bayes_R2(fit9_bayes, summary = T)
r2_bayes(fit9_bayes)
loo_fit8 <- loo(fit8_bayes)
loo_fit9 <- loo(fit9_bayes)
loo_compare(loo_fit8, loo_fit9)
posterior <- as.array(fit9_bayes)
#mcmc_pairs(posterior, pars = c('b_scalearea_p', 'b_scalemasa'))
mcmc_trace(posterior, pars = "b_scalesqrt_pl") + 
  xlab("Post-warmup iteration") 
mcmc_trace(posterior, pars = "b_scalewind.speed") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "b_scalebio12") + 
  xlab("Post-warmup iteration")
mcmc_trace(posterior, pars = "Intercept") + 
  xlab("Post-warmup iteration")
mcmc_acf(posterior, pars = "b_scalesqrt_pl", lags = 10)
mcmc_acf(posterior, pars = "b_scalewind.speed", lags = 10)
mcmc_acf(posterior, pars = "b_scalebio12", lags = 10)
mcmc_acf(posterior, pars = "Intercept", lags = 10)
mcmc_acf(posterior, pars = "shape", lags = 10)

mcmc_dens_overlay(posterior, pars = c("b_scalesqrt_pl", 'b_scalewind.speed',
                                      'b_scalebio12', 'Intercept','shape'))

pp = posterior_predict(fit9_bayes, ndraws = 1000)

y = df_sum2$sv
ppc_dens_overlay(y, pp)
pp_means <- apply(pp, 1, mean)  # Media de cada simulación
pp_sd <- apply(pp, 1, sd)       # Desviación estándar
pp_max <- apply(pp, 1, max)     # Máximo
pp_min <- apply(pp, 1, min)     # Mínimo

observed_mean <- mean(df_sum2$sv)
observed_sd <- sd(df_sum2$sv)
observed_max <- max(df_sum2$sv)
observed_min <- min(df_sum2$sv)



p_mean = ggplot(data.frame(pp_means), aes(x = pp_means)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_mean), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución de la Media Predicha", x = "Media Simulada de sv", y = "Frecuencia") + theme_minimal()


# Graficar la desviación estándar
p_sd = ggplot(data.frame(pp_sd), aes(x = pp_sd)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_sd), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución de la Desviación Estándar Predicha", x = "Desviación Estándar Simulada de sv", y = "Frecuencia")

# Graficar el máximo
p_max = ggplot(data.frame(pp_max), aes(x = pp_max)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_max), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución del Máximo Predicho", x = "Máximo Simulado de sv", y = "Frecuencia")

# Graficar el mínimo
p_min = ggplot(data.frame(pp_min), aes(x = pp_min)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  geom_vline(aes(xintercept = observed_min), color = "darkblue", linetype = "solid", linewidth = 1) +
  labs(title = "Distribución del Mínimo Predicho", x = "Mínimo Simulado de sv", y = "Frecuencia")

(p_mean+p_sd)/(p_max+p_min)



qres.fit9=createDHARMa(simulatedResponse = t(pp),
                       observedResponse = df_sum2$sv,
                       fittedPredictedResponse=apply(pp, 2, median))

s= data.frame(res = qnorm(residuals(qres.fit9)))

res.m9.brms=cbind(s,
                  df_sum2[,c("sqrt_pl", "wind.speed", "bio12")],
                  fitted=fitted(fit9_bayes, ndraws=1000)[,1],
                  pareto=loo(fit9_bayes, pointwise=T)$diagnostics$pareto_k)
p1 <- ggplot(res.m9.brms, aes(x = fitted, y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Fitted", y = "Residuals") +
  theme_minimal()

p2 <- ggplot(res.m9.brms, aes(x = scale(sqrt_pl), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "sqrt pl (std)", y = "Residuals") +
  theme_minimal()

p3 <- ggplot(res.m9.brms, aes(x = scale(wind.speed), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "wind speed (std)", y = "Residuals") +
  theme_minimal()
p4 <- ggplot(res.m9.brms, aes(x = scale(bio12), y = res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "temperatura (std)", y = "Residuals") +
  theme_minimal()
p5 <- ggplot(res.m9.brms, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()
p6 <- ggplot(res.m9.brms, aes(x = seq_along(pareto), y = pareto)) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Data point", y = "Pareto k") +
  theme_minimal()
library(patchwork)
(p1+p2+p3)/(p4+p5+p6)

### pareto k ######
loo_fit9 <- loo(fit9_bayes)
pareto_k <- loo_fit9$diagnostics$pareto_k
problematic_points <- which(pareto_k > 0.70)
spp = df_sum2[problematic_points,]$species
df_pl = df_sum2[!df_sum2$species%in%spp, ]
chk = name.check(tree2,df_pl)
tree_pl = drop.tip(tree2,chk$tree_not_data)
name.check(tree_pl,df_pl)

#binary_tree <- multi2di(tree2)
A <- ape::vcv.phylo(tree_pl, corr = T)
fit9_bayes_loo <- brm(
  sv ~  scale(sqrt_pl) + scale(wind.speed) + scale(bio12)  + (1 | gr(species, cov = A)),
  data = df_pl,
  data2 = list(A = A),  # Pass the covariance matrix here
  family = Gamma(link = "log"),
  prior = prior.m9,
  warmup=15000,chains=3, iter=30000, thin=5,
  control = list(max_treedepth = 25, stepsize = 0.001, adapt_delta = 0.999),
  cores = 11
)
summary(fit9_bayes_loo)
r2_bayes(fit9_bayes_loo)

## graficos final #####
preds = ggpredict(fit9_bayes,type="fixed", terms =
                    c("sqrt_pl[all]"))
preds_loo=ggpredict(fit9_bayes_loo,type="fixed", terms =
                      c("sqrt_pl[all]"))
p = ggplot(preds_loo, aes(x = x, y = predicted)) +
  geom_line(color = "grey35", linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = "grey35") +  # Banda de intervalo de confianza
  labs(x = "raíz cuadrada de PL", y = "Velocidad de asentamiento (m/s)") 
p_pl = p  + geom_line(data = preds, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = viridis(1)) + 
  geom_point(data = df_sum2, aes(x = sqrt_pl, y = sv), alpha = 0.2) + theme_bw()

preds_wind = ggpredict(fit9_bayes,type="fixed", terms =
                         c("wind.speed[all]"))
preds_wind_loo=ggpredict(fit9_bayes_loo,type="fixed", terms =
                           c("wind.speed[all]"))
p = ggplot(preds_wind_loo, aes(x = x, y = predicted)) +
  geom_line(color = "grey35", linetype = 'dashed') +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = "grey35") +  # Banda de intervalo de confianza
  labs(x = "Velocidad media del viento (m/s)", y = "Velocidad de asentamiento (m/s)") 
p_w = p  + geom_line(data = preds_wind, aes(x = x, y = predicted), color = viridis(1)) + 
  geom_ribbon(data = preds_wind, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = viridis(1)) + 
  geom_point(data = df_sum2, aes(x = wind.speed, y = sv), alpha = 0.2) + theme_bw() +  theme(axis.title.y =  element_blank())

preds_t_loo=ggpredict(fit9_bayes_loo,type="fixed", terms =
                        c("bio12[all]"))
preds_t_loo$model = 'modelo pareto K'
p_t = ggplot(preds_t_loo, aes(x = x, y = predicted)) +
  geom_line(aes(color = model, linetype = model)) +  # Línea de valores predichos
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) +  # Banda de intervalo de confianza
  labs(x = "Precipitación media anual (mm)", y = "Velocidad de asentamiento (m/s)") 

preds_t=ggpredict(fit9_bayes,type="fixed", terms =
                    c("bio12[all]"))
preds_t$model = 'modelo completo'
p_t = p_t + geom_line(data = preds_t, aes(x = x, y = predicted, color = model,linetype = model)) + 
  geom_ribbon(data = preds_t, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) + 
  geom_point(data = df_sum2, aes(x = bio12, y = sv), alpha = 0.2) +
  scale_color_manual(values = c(viridis(1),'grey35')) + 
  scale_fill_manual(values = c(viridis(1),'grey35')) + 
  scale_linetype_manual(values = c('solid','dashed'))+
  theme_bw() + 
  theme(axis.title.y =  element_blank(),
        legend.title = element_blank(), legend.position = c(1,1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.justification = c(1,1))
p_t
color_scheme_set('purple')
pi = mcmc_plot(fit9_bayes, prob=0.5, prob_outer = 0.95, regex_pars = c("^b")) + 
  geom_vline(aes(xintercept  = 0.0), linetype = 'solid', color = viridis(1, begin = 0.5)) + 
  scale_y_discrete(labels = c('intercepto', 'rc(PL)', 'vel viento', 'precipitación MA')) + 
  theme_bw() + 
  theme(axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())




combined_plot <- wrap_plots(wrap_plots(
  p_pl, p_w, p_t,ncol =3),               # Primera fila con 3 gráficos individuales
  wrap_plots(plot_spacer(),pi, plot_spacer(), ncol=3),  # Segunda fila combinada en 1 gráfico que ocupa toda la fila
  heights = c(2.5, 1)    # Distribución de alturas: primera fila más alta
) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ")")

pdf('modelo completo_pl.pdf',width = 10, height = 8, 'cm')
combined_plot
dev.off()

png('modelo completo_pl.png',width = 25, height = 15, units = 'cm', res = 600)

combined_plot
dev.off()
