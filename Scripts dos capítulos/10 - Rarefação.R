# Pacotes
library(iNEXT)
library(ecodados)
library(ggplot2)
library(vegan)
library(nlme)
library(dplyr)
library(piecewiseSEM)

# Dados
data("mite")
data("mite.xy")
coord <- mite.xy
colnames(coord) <- c("long", "lat")
data("mite.env")
agua <- mite.env[, 2]
dados_rarefacao <- ecodados::rarefacao_morcegos
rarefacao_repteis <- ecodados::rarefacao_repteis
rarefacao_anuros <- ecodados::rarefacao_anuros
dados_amostras <- ecodados::morcegos_rarefacao_amostras
 
# Curva de rarefação baseada no indivíduo (individual-based)
## Exemplo prático 1 - Morcegos

### Número de indivíduos por local
colSums(dados_rarefacao)

### Rarefação
# Datatype refere-se ao tipo de dados que você vai analisar (e.g. abundância, incidência).
# Endpoint refere-se ao valor máximo que você determina para a extrapolação.
resultados_morcegos <- iNEXT(dados_rarefacao, q = 0, 
                             datatype = "abundance", endpoint = 800)

### Gráfico
# type define o tipo de curva de rarefação
# 1 = curva de rarefação baseada no indivíduo ou amostra
# 2 = curva de representatividade da amostra
# 3 = curva de rarefação baseada na representatividade (coverage-based)

ggiNEXT(resultados_morcegos, type = 1) +
    geom_vline(xintercept = 166, lty = 2) +
    scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
    scale_colour_manual(values = c("darkorange", "darkorchid", "cyan4")) +
    scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
    labs(x = "Número de indivíduos", y = " Riqueza de espécies") +
    tema_livro()

## Exemplo prático 2 - Anuros e Répteis
### Análise
resultados_repteis <- iNEXT(rarefacao_repteis, q = 0,
                            datatype = "abundance", 
                            endpoint = 200)

### Visualizar os resultados
ggiNEXT(resultados_repteis, type = 1) +
  geom_vline(xintercept = 48, lty = 2) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  labs(x = "Número de indivíduos", y = " Riqueza de espécies") +
  tema_livro()

### Análise
resultados_anuros <- iNEXT(rarefacao_anuros, q = 0, 
                           datatype = "abundance", endpoint = 800)

### Visualizar os resultados
ggiNEXT(resultados_anuros, type = 1) + 
  geom_vline(xintercept = 37, lty = 2) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  scale_colour_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  labs(x = "Número de indivíduos", y = " Riqueza de espécies") +
  tema_livro()

### Curva de rarefação baseada em amostras (sample-based)
### Dados
# Usamos [,] para excluir os NAs. Lembrando que valores antes da 
# vírgula representam as linhas e os posteriores representam as colunas.
lista_rarefacao <- list(Tenentes = dados_amostras[1:18, 1],
                        Talhadinho = dados_amostras[, 2],
                        Experimental = dados_amostras[1:16, 3])

### Análise
res_rarefacao_amostras <- iNEXT(lista_rarefacao, q = 0, 
                                datatype = "incidence_freq")

### Gráfico
ggiNEXT(res_rarefacao_amostras , type = 1, color.var = "Assemblage") + 
  geom_vline(xintercept = 12, lty = 2) +
  scale_linetype_discrete(name = "Método", labels = c("Interpolado", "Extrapolado")) +
  scale_colour_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  labs(x = "Número de amostras", y = " Riqueza de espécies") +
  tema_livro()

# Curva de rarefação coverage-based
## Exemplo prático 4 - Morcegos
### Gráfico
# Visualizar os resultados da rarefação *coverage-based*. 
ggiNEXT(res_rarefacao_amostras, type = 3, color.var = "Assemblage") + 
  geom_vline(xintercept = 0.937, lty = 2) +
  scale_linetype_discrete(labels = c("Interpolado", "Extrapolado")) +
  scale_colour_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  labs(x = "Representatividade nas amostras", y = "Riqueza de espécies") +
  tema_livro()

## Exemplo prático 5 - Generalized Least Squares (GLS)
### Menos abundância
### Os dados estão com as comunidades nas colunas e as espécies nas linhas. 
### Para as análises teremos que transpor a planilha.
composicao_acaros <- as.data.frame(t(mite))

### Verificar qual é a menor abundância registrada nas comunidades. 
abun_min <- min(colSums(composicao_acaros))

### Riqueza rarefeita
resultados_rarefacao <- iNEXT(composicao_acaros, 
                              q = 0, 
                              datatype = "abundance", 
                              knots = abun_min,
                              endpoint = abun_min)
## Riqueza rarefeita
resultados_rarefacao <- iNEXT(composicao_acaros, 
                              q = 0, 
                              datatype = "abundance", 
                              knots = abun_min,
                              endpoint = abun_min)
### Riqueza rarefeita para cada comunidade (ATUALIZADO)
# Os resultados agora ficam agrupados no data.frame 'size_based'
resultados_comunidades <- resultados_rarefacao$iNextEst$size_based
# Filtramos a tabela para o tamanho alvo (m == abun_min) e pegamos os valores de qD (riqueza)
subset_res <- resultados_comunidades %>% 
  dplyr::filter(m == abun_min)
# O objeto subset_res já contém a riqueza rarefeita de TODAS as comunidades alinhadas!
riqueza_rarefeita <- subset_res$qD

### Dados finais
### Agrupando os dados em um data frame final.
dados_combinado <- data.frame(riqueza_rarefeita, agua, coord)

### Criando diferentes modelos usando a função gls
#### Sem estrutura espacial
no_spat_gls <- gls(riqueza_rarefeita ~ agua, data = dados_combinado, 
                   method = "REML")

#### Covariância esférica 
espher_model <- gls(riqueza_rarefeita ~ agua, data = dados_combinado, 
                    corSpher(form = ~lat + long, nugget = TRUE))

#### Covariância exponencial 
expon_model <- gls(riqueza_rarefeita ~ agua, data = dados_combinado, 
                   corExp(form = ~lat + long, nugget = TRUE))

#### Covariância Gaussiana 
gauss_model <- gls(riqueza_rarefeita ~ agua, data = dados_combinado, 
                   corGaus(form = ~lat + long, nugget = TRUE))

#### Covariância razão quadrática 
ratio_model <- gls(riqueza_rarefeita ~ agua, data = dados_combinado, 
                   corRatio(form = ~lat + long, nugget = TRUE))

### Seleção dos modelos
aic_fit <- AIC(no_spat_gls, espher_model, expon_model, 
               gauss_model, ratio_model)
aic_fit %>% arrange(AIC)

###
## Visualizando os resíduos do modelo selecionado
plot(gauss_model)

## Visualizando os resultados
summary(gauss_model)$tTable 


## Calculando o R-squared
rsquared(gauss_model)


## Obtendo os valores preditos pelo modelo
predito <- predict(gauss_model) 

## Plotando os resultados no gráfico
ggplot(data = dados_combinado, aes(x= agua, y= riqueza_rarefeita)) + 
  geom_point(size = 4, shape = 21, fill = "gray", alpha = 0.7) +
  geom_line(aes(y = predito), linewidth = 1) +
  labs(x = "Concentração de água no substrato", 
       y = "Riqueza rarefeita \ndas espécies de ácaros") +
  tema_livro()

