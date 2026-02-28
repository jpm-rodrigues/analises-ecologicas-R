# Pacotes
library(ecodados)
library(visdat)
library(tidyverse)
library(lattice)
library(RVAideMemoire)
library(DHARMa)
library(performance)
library(MuMIn)
library(piecewiseSEM)
library(MASS)
library(ggExtra)
library(Rmisc)
library(emmeans) 
library(sjPlot)
library(bbmle)
library(glmmTMB)
library(ordinal)
library(car)
library(ecolottery)
library(naniar)
library(vcd)
library(generalhoslem)

# Dados
lagartos <- ecodados::lagartos
parasitas <- ecodados::parasitas
fish <- ecodados::fish
fragmentos <- ecodados::fragmentos
uv_cells <- ecodados::uv_cells

# Dados de contagem: a distribuição de Poisson
## Gráfico
ggplot(fragmentos, aes(dfrag, Riqueza_obs)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm") +
  labs(x = "Distância para o fragmento mais próximo", 
       y = "Riqueza observada") +
  tema_livro()

## Modelo
mod_pois <- glm(Riqueza_obs ~ dfrag, family = poisson(link = "log"), data = fragmentos)

## Diagnose básica
par(mfrow = c(2, 2))
plot(mod_pois) 

## Diagnose avançada
simulationOutput <- simulateResiduals(fittedModel = mod_pois, plot = TRUE)

## Overdispersion
par(mfrow = c(1, 1))
testDispersion(mod_pois) # modelo tem overdispersion

## Testar a presença de overdispersion
check_overdispersion(mod_pois) # modelo tem overdispersion

## Resumo do modelo
summary(mod_pois)

## Dispersion parameter
deviance(mod_pois) / df.residual(mod_pois)


## Inflação de zeros - performanace
check_zeroinflation(mod_pois) # para diagnosticar se o modelo sofre de zero inflation


## Inflação de zeros - DHARMa
testZeroInflation(mod_pois) # para testar se existe zero inflation

## Coeficientes estimados pelo modelo
summary(mod_pois)

## Calculando o R2 do modelo
r.squaredGLMM(mod_pois)
r2(mod_pois)

## Plot do modelo predito
a1 <- ggplot(fragmentos, aes(dfrag, Riqueza_obs)) +
  geom_point(cex = 4,alpha = 0.7) +
  geom_smooth(method = "glm", formula = y~x, 
              method.args = list(family = "poisson"), se = TRUE) +
  labs(x = "Distância para o fragmento mais próximo", 
       y = "Riqueza observada") +
  tema_livro()

ggMarginal(a1, fill = "red")


## Ajuste do modelo
mod_nb <- glm.nb(Riqueza_obs ~ dfrag, data = fragmentos)

## Diagnose
par(mfrow = c(2, 2))
plot(mod_nb)
par(mfrow = c(1, 1))
(chat <- deviance(mod_nb) / df.residual(mod_nb)) # DISPERSION PARAMETER

## Diagnose avançada
simulationOutput <- simulateResiduals(fittedModel = mod_nb, plot = TRUE)

## Coeficiente de determinação
rsquared(mod_nb)

## Coeficientes estimados pelo modelo
summary(mod_nb)

## Plot do modelo
## Gráfico
ggplot(fragmentos, aes(dfrag, Riqueza_obs)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "glm.nb", formula = y~x, se = TRUE) +
  labs(x = "Distância para o fragmento mais próximo", 
       y = "Riqueza observada") +
  tema_livro()


# Dados de contagem: modelos quasi-likelihood
## Modelo
mod_quasipois <- glm(Riqueza_obs ~ dfrag, family = quasipoisson(link = "log"), data = fragmentos)

## Diagnose dos resíduos
EP <- resid(mod_quasipois, type = "pearson")
ED <- resid(mod_quasipois, type = "deviance")
mu <- predict(mod_quasipois, type = "response")
E <- fragmentos$Riqueza_obs - mu
EP2 <- E / sqrt(1.65662 * mu) # dispersion parameter da quasipoisson
op <- par(mfrow = c(2, 2))
plot(x = mu, y = E, main = "Response residuals")
plot(x = mu, y = EP, main = "Pearson residuals")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
plot(x = mu, y = ED, main = "Deviance residuals")
par(op)
par(mfrow = c(1, 1))

# Dados de contagem: a distribuição binomial
## Traduzir nomes das colunas e níveis de pigmentação 
colnames(uv_cells) <- c("UV", "Pigmentacao", "n_celulas", "linfocito", "neutrofilo", "basofilo", "monocito", "eosinofilo")
uv_cells$Pigmentacao[uv_cells$Pigmentacao=="Yes"] <- "sim"
uv_cells$Pigmentacao[uv_cells$Pigmentacao=="No"] <- "nao"

## Gráfico

# Calcular média e intervalo de confiança
eosinofilo <- summarySE(uv_cells, 
                        measurevar = "eosinofilo",
                        groupvars = c("UV", "Pigmentacao"))

# Definir posição de linhas e pontos no gráfico
pd <- position_dodge(0.1)

eosinofilo %>% 
  ggplot(aes(x = UV, y = eosinofilo, colour = Pigmentacao,
             group = Pigmentacao, fill = Pigmentacao)) +
  geom_errorbar(aes(ymin=eosinofilo-se, ymax=eosinofilo +se), 
                width=.1, size = 1.1, position=pd) +
  geom_line(position=pd, size = 1.1) +
  geom_point(pch = 21, colour = "black", position=pd, size=3.5) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  labs(x = "UV", y = "Eosinófilo", fill="Pigmentação", colour="Pigmentação")+
  tema_livro()

## Modelo
mod1 <- glm(cbind(eosinofilo, n_celulas) ~ UV * Pigmentacao, family = binomial, data = uv_cells)

## Diagnose dos resíduos
par(mfrow = c(2, 2))
plot(mod1)
par(mfrow = c(1, 1))

## Diagnose avançada
simulationBion <- simulateResiduals(fittedModel = mod1, plot = TRUE)

## Diagnose avançada
binned_residuals(mod1)

## Coeficientes estimados pelo modelo
summary(mod1)
anova(mod1)

## Parâmetros
pairs(emmeans(mod1, ~ UV|Pigmentacao))

ggplot(uv_cells, aes(UV, eosinofilo)) +
  geom_violin(aes(color = Pigmentacao)) +
  geom_jitter(shape = 16, position = position_jitter(0.1), cex = 4, alpha = 0.7) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  tema_livro()

# Análise com dados de incidência
## Traduzir nomes das colunas e níveis de pigmentação 
colnames(lagartos) <- c("numero", "sexo", "SVL", "comprimento_cauda", "cauda_autotomizada", "estado_cauda")
vis_miss(lagartos, cluster = TRUE) 

## Removendo dados faltantes
dados_semNA <- remove_missing(lagartos, vars = "sexo") 

## Visualizar
vis_miss(dados_semNA)

## Gráfico
ggplot(dados_semNA, aes(SVL, estado_cauda)) +
  geom_point(aes(shape = sexo, colour = sexo), size = 4, alpha = 0.4) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(y = "Estado da Cauda", x = "Comprimento Rostro-Cloacal (mm)", shape = "Sexo", colour = "Sexo") +
  tema_livro()

## Modelos
mod_log <- glm(estado_cauda ~ SVL * sexo, data = dados_semNA, family = binomial(link = "logit"))
mod_pro <- glm(estado_cauda ~ SVL * sexo, data = dados_semNA, family = binomial(link = "probit"))

# Seleção de modelos
AICctab(mod_log, mod_pro, nobs = 139)

## Diagnóse avançada
simulationBion <- simulateResiduals(fittedModel = mod_log, plot = T)

## Diagnóse avançada
binned_residuals(mod_log)

## Coeficientes estimados pelo modelo
summary(mod_log)
anova(mod_log, test = "Chisq" )

# Dados de contagem com excesso de zeros
## Explorando os dados com gráficos
ggplot(parasitas, aes(Raillietiella_mottae, fill = Especie)) +
  geom_density(alpha = 0.4) +
  facet_grid(Especie ~ Sexo) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  tema_livro() +
  theme(legend.position = "none")

ggplot(parasitas, aes(CRC, Raillietiella_mottae, fill = Especie)) +
  geom_point(size = 4, alpha = 0.4, shape = 21) +
  facet_grid(Sexo ~ Especie) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  theme(legend.position = "none") +
  labs(x = "Comprimento Rostro-Cloacal", y = expression(italic("Raillietiella mottae")))+
  tema_livro()

## Modelo
pois_plain <- glm(Raillietiella_mottae ~ CRC + Sexo * Especie, data = parasitas, family = "poisson")

## Diagnose avançada
# Verificar zero inflation
check_zeroinflation(pois_plain) # para diagnosticar se o modelo sofre de zero inflation
check_overdispersion(pois_plain)

## Modelos
# Hurdle model
hur_NB <- glmmTMB(Raillietiella_mottae ~ CRC + Sexo * Especie, zi = ~., data = parasitas, family = truncated_nbinom2) 

# zero-inflated Poisson
ziNB_mod2 <- glmmTMB(Raillietiella_mottae ~ CRC + Sexo * Especie, zi = ~., data = parasitas, family = nbinom2) 

# zero-inflated Negative binomial
ziP_mod2 <- glmmTMB(Raillietiella_mottae ~ CRC + Sexo * Especie, zi = ~., data = parasitas, family = poisson) 

## Diagnose de inflação de zeros
check_zeroinflation(hur_NB) 
check_zeroinflation(ziP_mod2)
check_zeroinflation(ziNB_mod2)

## Seleção de modelos
ICtab(pois_plain, hur_NB, ziP_mod2, ziNB_mod2, type = c("AICc"), weights = TRUE)

## Diagnoses Modelo Hurdle
simulationOutput <- simulateResiduals(fittedModel = hur_NB, plot = T)

## Diagnoses Modelo zero-inflated Poisson
simulationOutput <- simulateResiduals(fittedModel = ziP_mod2, plot = T)

## Diagnoses Modelo zero-inflated negative binomial
simulationOutput <- simulateResiduals(fittedModel = ziNB_mod2, plot = T) 

## Coeficientes estimados pelo modelo
summary(hur_NB)

## Gráfico
parasitas$phat <- predict(hur_NB, type = "response")
parasitas <- parasitas[with(parasitas, order(Sexo, Especie)), ]

ggplot(parasitas, aes(x = CRC, y = phat, colour = Especie,
                      shape = Sexo, linetype = Sexo)) +
  geom_point(aes(y = Raillietiella_mottae), size = 4, 
             alpha = .7, position = position_jitter(h = .2)) +
  geom_line(size = 1) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  labs(x = "Comprimento Rostro-Cloacal", 
       y = expression(paste("Abundância de ", italic("Raillietiella mottae")))) +
  tema_livro()

# Dados ordinais
cores <- read.csv2("https://ndownloader.figshare.com/files/10250700", header = TRUE)

## Tradução dos nomes das colunas
colnames(cores) <- c("animal", "tratamento", "tempo", "sexo", "preto", "vermelho")

## Filtrando dados - macho vermelho
macho_verm <- filter(cores, sexo == "M")

## Codificando como um fator ordenado 
macho_verm$animal <- factor(macho_verm$animal)
macho_verm$vermelho_ord <- factor(macho_verm$vermelho, 
                                  levels = c("1", "2", "3", "4", "5"), 
                                  ordered = TRUE)
str(macho_verm)

## Modelo
mod3 <- clmm(vermelho_ord ~ tratamento + tempo + (1|animal), data = macho_verm, threshold = "equidistant")

## Diagnose
assumption3 <- clm(vermelho_ord ~ tratamento + tempo, data = macho_verm, threshold = "equidistant")

scale_test(assumption3)
nominal_test(assumption3)

## Inferência
## Coeficientes estimados pelo modelo
summary(mod3)
anova(assumption3)
pairs(emmeans(mod3, ~ tratamento|tempo, adjust = "tukey"))

## Gráfico

# Calcular média e erro padrão
macho_verm_res <- summarySE(macho_verm, 
                            measurevar = "vermelho",
                            groupvars = c("tempo", "tratamento"))

# Definir posição de linhas e pontos no gráfico
pd <- position_dodge(0.1)

macho_verm_res %>% 
  ggplot(aes(x = tempo, y = vermelho, colour = tratamento,
             group = tratamento, fill = tratamento)) +
  geom_errorbar(aes(ymin=vermelho-se, ymax=vermelho +se), 
                width=.1, size = 1.1, position=pd) +
  geom_line(position=pd, size = 1.1) +
  geom_point(pch = 21, colour = "black", position=pd, size=3.5) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  xlab("Tempo de exposição (horas)") +
  ylab("Índice de eritróforos") +
  tema_livro()

# Dados contínuos: distribuição beta
## Tradução dos nomes das colunas
colnames(fish) <- c("animal", "tratamento", "tempo", "sexo", "preto", "vermelho")

## Filtrando os dados
fish$animal <- factor(fish$animal)
fish$sexo <- factor(fish$sexo)
macho_preto <- dplyr::filter(fish, sexo == "M")

## Gráfico
ggplot(macho_preto, aes(preto/100)) +
  geom_density(colour = "cyan4", fill = "cyan4", alpha = 0.4) +
  theme(legend.position = "none") +
  labs(x = "Índice de escuridão do corpo")+
  tema_livro()

## Modelo
mod2 <- glmmTMB(preto/100 ~ tratamento * tempo + (1|animal), family = beta_family, data = macho_preto)

## Diagnóse
simulationOutput <- simulateResiduals(fittedModel = mod2, plot = TRUE)

## Coeficientes estimados pelo modelo
Anova(mod2)

## níveis do fator da combinação
pairs(emmeans(mod2, ~ tratamento|tempo))

## Gráfico
escuridao <- summarySE(macho_preto, 
                       measurevar = "preto",
                       groupvars = c("tempo", "tratamento"))

# Definir posição de linhas e pontos no gráfico
pd <- position_dodge(0.1)

escuridao %>% 
  ggplot(aes(x = tempo, y = preto, colour = tratamento,
             group = tratamento, fill = tratamento)) +
  geom_errorbar(aes(ymin=preto-se, ymax=preto +se), 
                width=.1, size = 1.1, position=pd) +
  geom_line(position=pd, size = 1.1) +
  geom_point(pch = 21, colour = "black", position=pd, size=3.5) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  xlab("Tempo de experimento (horas)") +
  ylab("Índice de escuridão do corpo") +
  tema_livro()

