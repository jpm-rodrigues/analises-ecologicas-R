# ========= Pacotes ====================================
library(ecodados)
library(car)
library(ggpubr)
library(ggforce)
library(lsmeans)
library(lmtest)
library(sjPlot)
library(nlme)
library(ape)
library(fields)
library(tidyverse)
library(vegan)
library(rdist)

# ========= Dados ====================================
CRC_PN_macho <- ecodados::teste_t_var_igual
CRC_LP_femea <- ecodados::teste_t_var_diferente
Pareado <- ecodados::teste_t_pareado
correlacao_arbustos <- ecodados::correlacao
dados_regressao <- ecodados::regressoes
dados_regressao_mul <- ecodados::regressoes
dados_anova_simples <- ecodados::anova_simples
dados_dois_fatores <- ecodados::anova_dois_fatores
dados_dois_fatores_interacao <- ecodados::anova_dois_fatores
dados_dois_fatores_interacao2 <- ecodados::anova_dois_fatores_interacao2
dados_bloco <- ecodados::anova_bloco
dados_ancova <- ecodados::ancova
data("mite")
data("mite.xy")
coords <- mite.xy
colnames(coords) <- c("long", "lat")
data("mite.env")


# ========= Scripts =================================

# 7.1 - Teste T para duas amostras independentes

## Teste T para duas amostras com variâncias iguais

### Teste de normalidade
residuos <- lm(CRC ~ Estacao, data = CRC_PN_macho)
qqPlot(residuos)

### Teste de Shapiro-Wilk
residuos_modelo <- residuals(residuos)
shapiro.test(residuos_modelo)

### Teste de homogeneidade de variância
leveneTest(CRC ~ as.factor(Estacao), data = CRC_PN_macho)

### Análise Teste T 
t.test(CRC ~ Estacao, data = CRC_PN_macho, var.equal = TRUE)


### Gráfico
ggplot(data = CRC_PN_macho, aes(x = Estacao, y = CRC, color = Estacao)) + 
  labs(x = "Estações", 
       y = expression(paste("CRC (mm) - ", italic("P. nattereri")))) +
  geom_boxplot(fill = c("darkorange", "cyan4"), color = "black", 
               outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.7) +
  scale_color_manual(values = c("black", "black")) +
  tema_livro() +
  theme(legend.position = "none")

## Teste T para duas amostras independentes com variâncias diferentes

### Teste de normalidade usando QQ-plot
residuos_LP <- lm(CRC ~ Estacao, data = CRC_LP_femea)
qqPlot(residuos_LP)

### Teste de Shapiro-Wilk
residuos_modelo_LP <- residuals(residuos_LP)
shapiro.test(residuos_modelo_LP)

### Teste de homogeneidade da variância
leveneTest(CRC ~ as.factor(Estacao), data = CRC_LP_femea)

### Teste T
t.test(CRC ~ Estacao, data = CRC_LP_femea, var.equal = FALSE)

### Gráfico
ggplot(data = CRC_LP_femea, aes(x = Estacao, y = CRC, color = Estacao)) + 
  geom_boxplot(fill = c("darkorange", "cyan4"), width = 0.5, 
               color = "black", outlier.shape = NA, alpha = 0.7) +
  geom_jitter(shape = 20, position = position_jitter(0.2), color = "black", cex = 5) +
  scale_color_manual(values = c("darkorange", "cyan4")) +
  labs(x = "Estações", 
       y = expression(paste("CRC (mm) - ", italic("L. podicipinus"))), size = 15) +
  tema_livro() +
  theme(legend.position = "none")

## Teste T para amostras pareadas

### Análise Teste T Pareado
t.test(Riqueza ~ Estado, paired = TRUE, data = Pareado)

### Gráfico
ggpaired(Pareado, x = "Estado", y = "Riqueza",
         color = "Estado", line.color = "gray", line.size = 0.8, 
         palette = c("darkorange", "cyan4"), width = 0.5, 
         point.size = 4, xlab = "Estado das localidades", 
         ylab = "Riqueza de Espécies") +
  expand_limits(y = c(0, 150)) +
  tema_livro() 

# Correlação de Pearson

## Correlação de Pearson
cor.test(correlacao_arbustos$Tamanho_raiz, correlacao_arbustos$Tamanho_tronco, method = "pearson")

## Alternativamente
cor.test(~ Tamanho_tronco + Tamanho_raiz, data = correlacao_arbustos, method = "pearson")

## Gráfico
ggplot(data = correlacao_arbustos, aes(x = Tamanho_raiz, y = Tamanho_tronco)) + 
  labs(x = "Tamanho da raiz (m)", y = "Altura do tronco (m)") +
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  geom_text(x = 14, y = 14, label = "r = 0.89, P < 0.001", 
            color = "black", size = 5) +
  geom_smooth(method = lm, se = FALSE, color = "black", linetype = "dashed") +
  tema_livro() +
  theme(legend.position = "none")

# Regressão linear simples
## regressão simples
modelo_regressao <- lm(CRC ~ Temperatura, data = dados_regressao)

## Verificar as premissas do teste
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(modelo_regressao)
dev.off() # volta a configuração dos gráficos para o formato padrão

## Resultados usando a função anova
anova(modelo_regressao)

## Resultados usando a função summary
summary(modelo_regressao)

## Gráfico
ggplot(data = dados_regressao, aes(x = Temperatura, y = CRC)) + 
  labs(x = "Temperatura média anual (°C)", 
       y = "Comprimento rostro-cloacal (mm)") +
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  tema_livro() +
  theme(legend.position = "none")

# Regressão linear múltipla
## Regressão múltipla
modelo_regressao_mul <- lm(CRC ~ Temperatura + Precipitacao,
                           data = dados_regressao_mul)

# Multicolinearidade
vif(modelo_regressao_mul)

## Normalidade e homogeneidade das variâncias
plot_grid(plot_model(modelo_regressao_mul , type = "diag"))


## Resultados usando a função summary
summary(modelo_regressao_mul)


## Criando os modelos aninhados
modelo_regressao_mul <- lm(CRC ~ Temperatura + Precipitacao, 
                           data = dados_regressao_mul)
modelo_regressao <- lm(CRC ~ Temperatura, data = dados_regressao_mul)

## Likelihood-ratio test (LRT)
lrtest(modelo_regressao_mul, modelo_regressao)


## Comparando com o modelo somente com o intercepto
### Criando um modelo sem variáveis, só o intercepto.
modelo_intercepto <- lm(CRC ~ 1, data = dados_regressao_mul)
lrtest(modelo_regressao, modelo_intercepto)

# 7.6 Análises de Variância (ANOVA)
## ANOVA de um fator
### Análise ANOVA de um fator
Modelo_anova <- aov(Crescimento ~ Tratamento, data = dados_anova_simples) 


### Normalidade (shapiro)
shapiro.test(residuals(Modelo_anova))

### Homogeneidade da variância
bartlett.test(Crescimento ~ Tratamento, data = dados_anova_simples)

### Resultados da anova
anova(Modelo_anova)

### Diferenças entre os tratamentos
#### Teste de Tuckey's honest significant difference
TukeyHSD(Modelo_anova)


### Gráfico

#### Reorganizando a ordem que os grupos irão aparecer no gráfico
dados_anova_simples$Tratamento <- factor(dados_anova_simples$Tratamento,
                                         levels = c("Controle", "Adubo_Tradicional", "Adubo_X-2020"))


ggplot(data = dados_anova_simples, 
       aes(x = Tratamento, y = Crescimento, color = Tratamento)) + 
  geom_boxplot(fill = c("darkorange", "darkorchid", "cyan4"), 
               color = "black", show.legend = FALSE, alpha = 0.4) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
  geom_text(x = 1, y = 12, label = "ab", color = "black", size = 5) +
  geom_text(x = 2, y = 17, label = "a", color = "black", size = 5) +
  geom_text(x = 3, y = 17, label = "b", color = "black", size = 5) +
  scale_x_discrete(labels = c("Sem adubo", "Tradicional", "X-2020")) +
  labs(x = "Adubação", y = "Crescimento Coffea arabica (cm)", size = 20) +
  tema_livro() +
  theme(legend.position = "none") 


## ANOVA com dois fatores

# A interação entre os fatores é representada por *
Modelo1 <- aov(Tempo ~ Pessoas * Idade, data = dados_dois_fatores) 

# Olhando os resultados
anova(Modelo1)

# Criando modelo sem interação.
Modelo2 <- aov(Tempo ~ Pessoas + Idade, data = dados_dois_fatores) 

# Verificando as premissas do teste.
plot_grid(plot_model(Modelo2, type = "diag"))

# Resultados do modelo
anova(Modelo2)

## Gráfico
ggplot(data = dados_dois_fatores_interacao, 
       aes(y = Tempo, x = Pessoas, color = Idade)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom ="point", aes(group = Idade, x = Pessoas),
               color = "black",
               position = position_dodge(0.7), size  = 4) +
  geom_link(aes(x = 0.8, y = 31, xend = 1.8, yend = 40), color = "darkorange", 
            lwd  = 1.3, linetype = 2) + 
  geom_link(aes(x = 1.2, y = 19, xend = 2.2, yend = 26.5), 
            color = "cyan4", lwd  = 1.3, linetype = 2) + 
  labs(x = "Sistema XY de determinação do sexo", 
       y = "Tempo (horas) para eliminar a droga") +
  scale_color_manual(values = c("darkorange", "cyan4", "darkorange", "cyan4")) +
  scale_y_continuous(limits = c(10, 50), breaks = c(10, 20, 30, 40, 50)) +
  tema_livro()  

## ANOVA com dois fatores com efeito da interação

### Análise anova de dois fatores 
Modelo_interacao2 <- aov(Tempo ~ Pessoas * Idade, 
                         data = dados_dois_fatores_interacao2)

### Olhando os resultados
anova(Modelo_interacao2)

### Gráfico
ggplot(data = dados_dois_fatores_interacao2, 
       aes(y = Tempo, x = Pessoas, color = Idade)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom ="point", aes(group = Idade, x = Pessoas), 
               color = "black", position = position_dodge(0.7), size  = 4) +
  geom_link(aes(x = 0.8, y = 31, xend = 1.8, yend = 27), color = "darkorange", 
            lwd  = 1.3, linetype = 2) + 
  geom_link(aes(x = 1.2, y = 19, xend = 2.2, yend = 41), color = "cyan4", 
            lwd  = 1.3, linetype = 2) + 
  labs(x = "Sistema XY de determinação do sexo", 
       y = "Tempo (horas) para eliminar a droga") +
  scale_color_manual(values = c("darkorange", "cyan4", "darkorange", "cyan4")) +
  scale_y_continuous(limits = c(10, 50), breaks = c(10, 20, 30, 40, 50)) +
  tema_livro() 


## ANOVA em blocos aleatorizados 

### Análise Anova em blocos aleatorizados
model_bloco <- aov(Riqueza ~ Pocas + Error(Blocos), data = dados_bloco)
summary(model_bloco)

##$ Forma errada de análisar Anova em blocos
modelo_errado <- aov(Riqueza ~ Pocas, data = dados_bloco)
anova(modelo_errado)

### Teste de Tuckey's honest significant difference
pairs(lsmeans(model_bloco, "Pocas"), adjust = "tukey")

### Gráfico

#### Reordenando a ordem que os grupos irão aparecer no gráfico. ####
dados_bloco$Pocas <- factor(dados_bloco$Pocas, 
                            levels = c("Int-100m", "Int-50m", "Borda", "Mat-50m", "Mat-100m"))

#### 

ggplot(data = dados_bloco, aes(x = Pocas, y = Riqueza)) + 
  labs(x = "Poças artificiais", y = "Riqueza de espécies de anuros") +
  geom_boxplot(color = "black", show.legend = FALSE, alpha = 0.4) +
  geom_jitter(shape = 16, position = position_jitter(0.1), cex = 4, alpha = 0.7) +
  scale_x_discrete(labels = c("-100m","-50m","Borda", "50m", "100m")) +
  tema_livro() +
  theme(legend.position = "none") 

### Análise de covariância (ANCOVA)


### Ancova
modelo_ancova <- lm(Biomassa ~ Herbivoria * Raiz, data = dados_ancova)

### Verificando as premissas da Ancova
plot_grid(plot_model(modelo_ancova, type = "diag"))

### Resultados do modelo
anova(modelo_ancova)

### Criando modelo sem interação
modelo_ancova2 <- lm(Biomassa ~ Herbivoria + Raiz, data = dados_ancova)

### Likelihood-ratio test
lrtest(modelo_ancova, modelo_ancova2)


### Gráfico
ggplot(data = dados_ancova, aes(x = Raiz, y = Biomassa, fill = Herbivoria)) + 
  labs(x = "Tamanho da raiz (cm)", y = "Biomassa dos frutos (g)") +
  geom_point(size = 4, shape = 21, alpha = 0.7) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4"),
                    labels = c("Com herbivoria", "Sem herbivoria")) +
  geom_smooth(aes(color = Herbivoria), method = "lm", show.legend = FALSE) +
  tema_livro()


# Generalized Least Squares (GLS)

### Calcular a riqueza de espécies em cada comunidade
riqueza <- specnumber(mite) 

### Selecionar a variável ambiental - quantidade de água no substrato
agua <- mite.env[,2]

### Criar um data.frame com riqueza, quantidade de água no substrato e coordenadas geográficas
mite_dat <- data.frame(riqueza, agua, coords)

### Modelo linear sem incorporar a estrutura espacial
#### Modelo
linear_model <- lm(riqueza ~ agua, mite_dat) 

#### Resíduos
par(mfrow = c(2, 2)) 
plot(linear_model, which = 1:4)

#### Resultados do modelo
res_lm <- summary(linear_model)

#### Coeficiente de determinação e coeficientes
res_lm$adj.r.squared
res_lm$coefficients

#### Modelo gls sem estrutura espacial
no_spat_gls <- gls(riqueza ~ agua, mite_dat, method = "REML")

#### Variograma
variog_mod1 <- nlme::Variogram(no_spat_gls, form = ~lat+long, 
                               resType = "normalized") 

#### Gráfico
dev.off() # reseta o gráfico
plot(variog_mod1)

#### Primeiro precisamos calcular uma matriz de distâncias geográficas entre as comunidades
dat_dist <- pdist(coords) # matriz de distância

Moran.I(x = mite_dat$riqueza, w = dat_dist)


#### Covariância esférica
espher_model <- gls(riqueza ~ agua, mite_dat, 
                    corSpher(form = ~lat+long, nugget = TRUE))

#### Covariância exponencial
expon_model <- gls(riqueza ~ agua, mite_dat, 
                   corExp(form = ~lat+long, nugget = TRUE))

#### Covariância Gaussiana
gauss_model <- gls(riqueza ~ agua, mite_dat, 
                   corGaus(form = ~lat+long, nugget = TRUE))

#### Covariância linear
cor_linear_model <- gls(riqueza ~ agua, mite_dat, 
                        corLin(form = ~lat+long, nugget = TRUE))

#### Covariância razão quadrática
ratio_model <- gls(riqueza ~ agua, mite_dat, 
                   corRatio(form = ~lat+long, nugget = TRUE))

#### Seleção de modelos
aic_fit <- AIC(no_spat_gls, espher_model, expon_model, gauss_model, cor_linear_model, ratio_model)
aic_fit %>% arrange(AIC)

#### Gráfico
dev.off()
plot(residuals(ratio_model, type = "normalized") ~ fitted(ratio_model))


#### Varigrama
ratio_variog <- Variogram(ratio_model, form = ~lat+long, resType = "normalized")

#### Resumo dos modelos
summary(ratio_model)$tTable
summary(no_spat_gls)$tTable

#### Gráficos
plot(ratio_variog, main = "Variograma como Modelo Ratio")
plot(variog_mod1, main = "Variograma Modelo Normal")
