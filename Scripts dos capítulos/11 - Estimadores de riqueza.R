## Pacotes
library(iNEXT)
library(devtools)
library(ecodados)
library(ggplot2)
library(vegan)
library(nlme)
library(dplyr)
library(piecewiseSEM)

## Dados
dados_coleta <- poca_anuros
data(mite)
data(mite.xy)
coord <- mite.xy
colnames(coord) <- c("long", "lat") # altera o nome das colunas
data(mite.env)
agua <- mite.env[, 2] # seleciona a variável de interesse

# CHAO 1 - Exemplo prático

## Análise
est_chao1 <- estaccumR(dados_coleta, permutations = 100)
summary(est_chao1, display = "chao")

## Preparando os dados para fazer o gráfico
resultados <- summary(est_chao1, display = c("S", "chao"))
res_chao <- cbind(resultados$chao[, 1:4], resultados$S[, 2:4])
res_chao <- as.data.frame(res_chao)
colnames(res_chao) <- c(
  "Amostras", "Chao", "C_inferior", "C_superior",
  "Riqueza", "R_inferior", "R_superior"
)

## Gráfico
ggplot(res_chao, aes(y = Riqueza, x = Amostras)) +
  geom_point(aes(y = Chao, x = Amostras + 0.1),
    size = 4,
    color = "darkorange", alpha = 0.7
  ) +
  geom_point(aes(y = Riqueza, x = Amostras),
    size = 4,
    color = "cyan4", alpha = 0.7
  ) +
  geom_point(y = 7.5, x = 9, size = 4, color = "darkorange", alpha = 0.7) +
  geom_point(y = 5.9, x = 9, size = 4, color = "cyan4", alpha = 0.7) +
  geom_label(y = 7.5, x = 12, label = "Riqueza estimada - Chao 1") +
  geom_label(y = 5.9, x = 11.3, label = "Riqueza observada") +
  geom_line(aes(y = Chao, x = Amostras), color = "darkorange") +
  geom_line(aes(y = Riqueza, x = Amostras), color = "cyan4") +
  geom_linerange(aes(
    ymin = C_inferior, ymax = C_superior,
    x = Amostras + 0.1
  ), color = "darkorange") +
  geom_linerange(aes(
    ymin = R_inferior, ymax = R_superior,
    x = Amostras
  ), color = "cyan4") +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 1)) +
  labs(x = "Número de amostras", y = "Riqueza estimada - Chao 1") +
  tema_livro()

# ACE - Abundance-based Coverage Estimator
## Análise
est_ace <- estaccumR(dados_coleta, permutations = 100)
summary(est_ace, display = "ace")

## Preparando os dados para fazer o gráfico
resultados_ace <- summary(est_ace, display = c("S", "ace"))
res_ace <- cbind(resultados_ace$ace[, 1:4], resultados_ace$S[, 2:4])
res_ace <- as.data.frame(res_ace)
colnames(res_ace) <- c(
  "Amostras", "ACE", "ACE_inferior", "ACE_superior",
  "Riqueza", "R_inferior", "R_superior"
)

## Gráfico
ggplot(res_ace, aes(y = Riqueza, x = Amostras)) +
  geom_point(aes(y = ACE, x = Amostras + 0.1),
    size = 4,
    color = "darkorange", alpha = 0.7
  ) +
  geom_point(aes(y = Riqueza, x = Amostras),
    size = 4,
    color = "cyan4", alpha = 0.7
  ) +
  geom_point(y = 7.5, x = 9, size = 4, color = "darkorange", alpha = 0.7) +
  geom_point(y = 5.9, x = 9, size = 4, color = "cyan4", alpha = 0.7) +
  geom_label(y = 7.5, x = 11.7, label = "Riqueza estimada - ACE") +
  geom_label(y = 5.9, x = 11.3, label = "Riqueza observada") +
  geom_line(aes(y = ACE, x = Amostras), color = "darkorange") +
  geom_line(aes(y = Riqueza, x = Amostras), color = "cyan4") +
  geom_linerange(aes(
    ymin = ACE_inferior, ymax = ACE_superior,
    x = Amostras + 0.1
  ), color = "darkorange") +
  geom_linerange(aes(
    ymin = R_inferior, ymax = R_superior,
    x = Amostras
  ), color = "cyan4") +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 1)) +
  labs(x = "Número de amostras", y = "Riqueza estimada - ACE") +
  tema_livro()

# Estimadores baseados na incidência das espécies
# CHAO 2 - Exemplo prático
## Análise
est_chao2 <- poolaccum(dados_coleta, permutations = 100)
summary(est_chao2, display = "chao")

## Preparando os dados para fazer o gráfico
resultados_chao2 <- summary(est_chao2, display = c("S", "chao"))
res_chao2 <- cbind(resultados_chao2$chao[, 1:4], resultados_chao2$S[, 2:4])
res_chao2 <- as.data.frame(res_chao2)
colnames(res_chao2) <- c(
  "Amostras", "Chao2", "C_inferior", "C_superior",
  "Riqueza", "R_inferior", "R_superior"
)

## Gráfico
ggplot(res_chao2, aes(y = Riqueza, x = Amostras)) +
  geom_point(aes(y = Chao2, x = Amostras + 0.1),
    size = 4,
    color = "darkorange", alpha = 0.7
  ) +
  geom_point(aes(y = Riqueza, x = Amostras),
    size = 4,
    color = "cyan4", alpha = 0.7
  ) +
  geom_point(y = 9.8, x = 10, size = 4, color = "darkorange", alpha = 0.7) +
  geom_point(y = 7.7, x = 10, size = 4, color = "cyan4", alpha = 0.7) +
  geom_label(y = 9.8, x = 12.95, label = "Riqueza estimada - Chao 2") +
  geom_label(y = 7.7, x = 12.3, label = "Riqueza observada") +
  geom_line(aes(y = Chao2, x = Amostras), color = "darkorange") +
  geom_line(aes(y = Riqueza, x = Amostras), color = "cyan4") +
  geom_linerange(aes(
    ymin = C_inferior, ymax = C_superior,
    x = Amostras + 0.1
  ), color = "darkorange") +
  geom_linerange(aes(
    ymin = R_inferior, ymax = R_superior,
    x = Amostras
  ), color = "cyan4") +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 1)) +
  labs(x = "Número de amostras", y = "Riqueza estimada - Chao 2") +
  tema_livro()

# JACKKNIFE 1
## Análise
est_jack1 <- poolaccum(dados_coleta, permutations = 100)
summary(est_jack1, display = "jack1")

## Preparando os dados para fazer o gráfico
resultados_jack1 <- summary(est_jack1, display = c("S", "jack1"))
res_jack1 <- cbind(resultados_jack1$jack1[, 1:4], resultados_jack1$S[, 2:4])
res_jack1 <- as.data.frame(res_jack1)
colnames(res_jack1) <- c(
  "Amostras", "JACK1", "JACK1_inferior", "JACK1_superior",
  "Riqueza", "R_inferior", "R_superior"
)

## Gráfico
ggplot(res_jack1, aes(y = Riqueza, x = Amostras)) +
  geom_point(aes(y = JACK1, x = Amostras + 0.1),
    size = 4,
    color = "darkorange", alpha = 0.7
  ) +
  geom_point(aes(y = Riqueza, x = Amostras),
    size = 4,
    color = "cyan4", alpha = 0.7
  ) +
  geom_point(y = 9.9, x = 9, size = 4, color = "darkorange", alpha = 0.7) +
  geom_point(y = 8.6, x = 9, size = 4, color = "cyan4", alpha = 0.7) +
  geom_label(y = 9.9, x = 12.5, label = "Riqueza estimada - Jackknife 1") +
  geom_label(y = 8.6, x = 11.5, label = "Riqueza observada") +
  geom_line(aes(y = JACK1, x = Amostras), color = "darkorange") +
  geom_line(aes(y = Riqueza, x = Amostras), color = "cyan4") +
  geom_linerange(aes(
    ymin = JACK1_inferior, ymax = JACK1_superior,
    x = Amostras + 0.1
  ), color = "darkorange") +
  geom_linerange(aes(
    ymin = R_inferior, ymax = R_superior,
    x = Amostras
  ), color = "cyan4") +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 1)) +
  labs(x = "Número de amostras", y = "Riqueza estimada - Jackknife 1") +
  tema_livro()

# JACKKNIFE 2
## Análise
est_jack2 <- poolaccum(dados_coleta, permutations = 100)
summary(est_jack2, display = "jack2")

## Preparando os dados para fazer o gráfico
resultados_jack2 <- summary(est_jack2, display = c("S", "jack2"))
res_jack2 <- cbind(resultados_jack2$jack2[, 1:4], resultados_jack2$S[, 2:4])
res_jack2 <- as.data.frame(res_jack2)
colnames(res_jack2) <- c(
  "Amostras", "JACK2", "JACK2_inferior", "JACK2_superior",
  "Riqueza", "R_inferior", "R_superior"
)

## Gráfico
ggplot(res_jack2, aes(y = Riqueza, x = Amostras)) +
  geom_point(aes(y = JACK2, x = Amostras + 0.1),
    size = 4,
    color = "darkorange", alpha = 0.7
  ) +
  geom_point(aes(y = Riqueza, x = Amostras),
    size = 4,
    color = "cyan4", alpha = 0.7
  ) +
  geom_point(y = 9.9, x = 9, size = 4, color = "darkorange", alpha = 0.7) +
  geom_point(y = 8.2, x = 9, size = 4, color = "cyan4", alpha = 0.7) +
  geom_label(y = 9.9, x = 12.5, label = "Riqueza estimada - Jackknife 2") +
  geom_label(y = 8.2, x = 11.5, label = "Riqueza observada") +
  geom_line(aes(y = JACK2, x = Amostras), color = "darkorange") +
  geom_line(aes(y = Riqueza, x = Amostras), color = "cyan4") +
  geom_linerange(aes(
    ymin = JACK2_inferior, ymax = JACK2_superior,
    x = Amostras + 0.1
  ), color = "darkorange") +
  geom_linerange(aes(
    ymin = R_inferior, ymax = R_superior,
    x = Amostras
  ), color = "cyan4") +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 1)) +
  labs(x = "Número de amostras", y = "Riqueza estimada - Jackknife 2") +
  tema_livro()

# BOOTSTRAP
est_boot <- poolaccum(dados_coleta, permutations = 100)
summary(est_boot, display = "boot")

## Preparando os dados para fazer o gráfico
resultados_boot <- summary(est_boot, display = c("S", "boot"))
res_boot <- cbind(resultados_boot$boot[, 1:4], resultados_boot$S[, 2:4])
res_boot <- as.data.frame(res_boot)
colnames(res_boot) <- c(
  "Amostras", "BOOT", "BOOT_inferior", "BOOT_superior",
  "Riqueza", "R_inferior", "R_superior"
)
## Gráfico
ggplot(res_boot, aes(y = Riqueza, x = Amostras)) +
  geom_point(aes(y = BOOT, x = Amostras + 0.1),
    size = 4,
    color = "darkorange", alpha = 0.7
  ) +
  geom_point(aes(y = Riqueza, x = Amostras),
    size = 4,
    color = "cyan4", alpha = 0.7
  ) +
  geom_point(y = 10.4, x = 9, size = 4, color = "darkorange", alpha = 0.7) +
  geom_point(y = 9.3, x = 9, size = 4, color = "cyan4", alpha = 0.7) +
  geom_label(y = 10.4, x = 12.3, label = "Riqueza estimada - Bootstrap") +
  geom_label(y = 9.3, x = 11.5, label = "Riqueza observada") +
  geom_line(aes(y = BOOT, x = Amostras), color = "darkorange") +
  geom_line(aes(y = Riqueza, x = Amostras), color = "cyan4") +
  geom_linerange(aes(
    ymin = BOOT_inferior, ymax = BOOT_superior,
    x = Amostras + 0.1
  ), color = "darkorange") +
  geom_linerange(aes(
    ymin = R_inferior, ymax = R_superior,
    x = Amostras
  ), color = "cyan4") +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 1)) +
  labs(x = "Número de amostras", y = "Riqueza estimada - Bootstrap") +
  tema_livro()

# Interpolação e extrapolação baseadas em rarefação usando amostragens de incidência ou abundância
# Exercício 1
## Preparando os dados para análises considerando a abundância
dados_inext_abu <- colSums(dados_coleta)
resultados_abundancia <- iNEXT(dados_inext_abu,
  q = 0, datatype = "abundance",
  endpoint = 600
)

## Gráfico
ggiNEXT(resultados_abundancia, type = 1) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  scale_colour_manual(values = "darkorange") +
  scale_fill_manual(values = "darkorange") +
  labs(x = "Número de indivíduos", y = " Riqueza de espécies") +
  tema_livro()

## Interpretando os resultados
### Preparando os dados para análises considerando a incidência
#### Precisa transpor o data frame.
dados_inext <- as.incfreq(t(dados_coleta))
resultados_incidencia <- iNEXT(dados_inext,
  q = 0, datatype = "incidence_freq",
  endpoint = 28
)

#### Gráfico
ggiNEXT(resultados_incidencia, type = 1, color.var = "Order.q") +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  scale_colour_manual(values = "darkorange") +
  scale_fill_manual(values = "darkorange") +
  labs(x = "Número de amostras", y = " Riqueza de espécies") +
  tema_livro()

## Exercício 2
## Transposição
# Os dados estão com as comunidades nas colunas e as espécies nas linhas.
# Para as análises teremos que transpor a planilha.
composicao_acaros <- as.data.frame(t(mite))

## Dados
# A comunidade com maior abundância tem 781 indivíduos.
max(colSums(composicao_acaros))
#> [1] 781

# Calcular a riqueza extrapolada de espécies para todas as comunidades
# considerando a maior abundância.
## Essa demora alguns minutos para rodar
resultados_extrapolacao <- iNEXT(composicao_acaros,
  q = 0,
  datatype = "abundance",
  endpoint = 781
)

# O iNEXT v3 salva os resultados consolidados de rarefação/extrapolação no objeto size_based
resultados_comunidades_ext <- resultados_extrapolacao$iNextEst$size_based

## Extraindo os Resultados
# O objetivo é extrair a riqueza estimada extrapolada para 781 individuos
# Utilizando as funções do dplyr (já carregado) para manter a ordem original das comunidades
riqueza_extrapolada <- resultados_comunidades_ext %>%
  filter(m == 781) %>%
  pull(qD)

## Dados
# Criando data frame com todas as variáveis
dados_combinado_ext <- data.frame(riqueza_extrapolada, agua, coord)

## Modelo gls sem estrutura espacial
no_spat_gls <- gls(riqueza_extrapolada ~ agua, 
                   data = dados_combinado_ext, 
                   method = "REML")

## Covariância esférica
espher_model <- gls(riqueza_extrapolada ~ agua, 
                    data = dados_combinado_ext, 
                    corSpher(form = ~lat + long, nugget = TRUE))

## Covariância exponencial
expon_model <- gls(riqueza_extrapolada ~ agua, 
                   data = dados_combinado_ext, 
                   corExp(form = ~lat + long, nugget = TRUE))

## Covariância Gaussiana
gauss_model <- gls(riqueza_extrapolada ~ agua,
                   data = dados_combinado_ext, 
                   corGaus(form = ~lat + long, nugget = TRUE))

## Covariância razão quadrática
ratio_model <- gls(riqueza_extrapolada ~ agua, 
                   data = dados_combinado_ext, 
                   corRatio(form = ~lat + long, nugget = TRUE))
## Seleção dos modelos
aic_fit_ext <- AIC(no_spat_gls, espher_model, expon_model, gauss_model, ratio_model)
aic_fit_ext %>% arrange(AIC)

## Visualizando os resíduos do modelo com menor valor de AIC (veja Capítulo 7)
plot(ratio_model)

## Visualizando os resultados e calculando pseudo-R-squared.
summary(ratio_model)$tTable 
#>                   Value   Std.Error   t-value      p-value
#> (Intercept) 24.09577588 4.816461582  5.002796 4.227862e-06
#> agua        -0.01181425 0.006977381 -1.693221 9.499017e-02
rsquared(ratio_model)
#>              Response   family     link method  R.squared
#> 1 riqueza_extrapolada gaussian identity   none 0.05977552

## Gráfico
predito_ext <- predict(ratio_model) 

ggplot(data = dados_combinado_ext, aes(x= agua, y= riqueza_extrapolada)) + 
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  geom_line(aes(y = predito_ext), size = 1) +
  labs(x = "Concentração de água no substrato", 
       y = "Riqueza extrapolada \ndas espécies de ácaros") +
  tema_livro()
