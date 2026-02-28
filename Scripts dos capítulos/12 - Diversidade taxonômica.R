## Pacotes
library(devtools)
library(ecodados)
library(vegan)
library(ggplot2)
library(BiodiversityR)
library(hillR)
library(betapart)

## Dados
composicao_especies <- ecodados::composicao_anuros_div_taxonomica
precipitacao <- ecodados::precipitacao_div_taxonomica

# Diversidade alfa

## Calculando a riqueza observada de espécies para cada comunidade
riqueza_sp <- specnumber(composicao_especies)
riqueza_sp

## Calculamos a abundância total para cada comunidade
abundancia <- apply(composicao_especies, 1, sum)
abundancia

## Índice de Margalef
# A função round é para limitar o resultado para duas casas decimais.
Margalef <- round((riqueza_sp - 1)/log(abundancia), 2)
Margalef

## Índice de Menhinick
Menhinick <- round(riqueza_sp/sqrt(abundancia), 2)
Menhinick

## Juntando todos os dados em um único data frame
dados <- data.frame(precipitacao$prec, riqueza_sp, Margalef, Menhinick)

## Renomenado as colunas
colnames(dados) <- c("Precipitacao", "Riqueza", "Margalef", "Menhinick")

## ANOVA (riqueza de espécies e a precipitação anual.)
anova_riq <- lm(Riqueza ~ Precipitacao, data = dados)
anova(anova_riq)


## ANOVA (Margalef e a precipitação anual.)
anova_marg <- lm(Margalef ~ Precipitacao, data = dados)
anova(anova_marg)

## ANOVA (Menhinick e a precipitação anual)
anova_menh <- lm(Menhinick ~ Precipitacao, data = dados)
anova(anova_menh)

## Gráfico
ggplot(data = dados, aes(x= Precipitacao, y= Riqueza)) + 
  labs(x = "Precipitação anual (mm)", y = "Riqueza de espécies") +
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  tema_livro()

## Diversidade de espécies

## Índice de Shannon
# MARGIN = 1 significa que a função irá calcular o índice considerando 
# as linhas do data.frame (comunidades).
shannon_res <- diversity(composicao_especies, index = "shannon", MARGIN = 1)
shannon_res

## Índice de Simpson
simpson_res <- diversity(composicao_especies, index = "simpson", MARGIN = 1) 
simpson_res

## Índice de Pielou
Pielou <- shannon_res/log(specnumber(composicao_especies))

## Juntando todos os dados em um único data frame
dados_div <- data.frame(precipitacao$prec, shannon_res, 
                        simpson_res, Pielou)

## Renomeando as colunas
colnames(dados_div) <- c("Precipitacao", "Shannon", "Simpson", "Pielou")
 
## ANOVA (shannon e precipitação)
anova_shan <- lm(Shannon ~ Precipitacao, data = dados_div)
anova(anova_shan)

## ANOVA (simpson e precipitação)
anova_simp <- lm(Simpson ~ Precipitacao, data = dados_div)
anova(anova_simp)

## ANOVA (pielou e precipitação)
anova_piel <- lm(Pielou ~ Precipitacao, data = dados_div)
anova(anova_piel)



# Diagramas de Whittaker ou Curva de Dominância
## Cálculo da curva para as comunidades 2 e 3
rank_com2 <- rankabundance(composicao_especies[2, composicao_especies[2,] > 0])
rank_com3 <- rankabundance(composicao_especies[3, composicao_especies[3,] > 0])

## Gráfico (diagrama de whittaker)
## Veja a ajuda da função rankabundplot para outros exemplos de gráficos.
rankabunplot(rank_com2, scale = "logabun", specnames = c(1), 
             pch = 19, col = "darkorange")
rankabunplot(rank_com3, scale = "logabun", specnames = c(1), pch = 19, 
             xlim = c(0,10), addit = TRUE, col = "cyan4")
legend(5, 40, legend = c("Comunidade 2", "Comunidade 3"),
       col = c("darkorange", "cyan4"), lty = 1, cex = 0.8, box.lty = 0)


# Curvas de distribuição de abundâncias
## Teste das curvas de distribuição de abundâncias - Comunidade 2
curvas_dominancia_com2 <- radfit(composicao_especies[2,])
curvas_dominancia_com2
plot(curvas_dominancia_com2, 
     ylab = "Abundância", 
     xlab = "Ranqueamento das espécies")
## Interpretação dos resultados 
## Teste das curvas de distribuição de abundâncias
curvas_dominancia_todas <- radfit(composicao_especies)
curvas_dominancia_todas
# Vamos fazer um gráfico para cada comunidade
plot(curvas_dominancia_todas, log = "y")

# Números de Hill ou Série de Hill
## Número de Hill para q = 0
hill_res_q_0 <- hill_taxa(composicao_especies, q  = 0)
hill_res_q_0

## Número de Hill para q = 1
hill_res_q_1 <- hill_taxa(composicao_especies, q  = 1)
hill_res_q_1

## Número de Hill para q = 2
hill_res_q_2 <- hill_taxa(composicao_especies, q  = 2)
hill_res_q_2

## Resultados
res_hill <- data.frame(hill_res_q_0, hill_res_q_1, hill_res_q_2)
colnames(res_hill) <- c("q=0", "q=1", "q=2")
head(res_hill)

# Diversidade beta
## Transformando dados em presencia e ausência.
composicao_PA <- decostand(composicao_especies, method = "pa")

## Diversidade beta
resultado_PA <- beta.pair(composicao_PA, index.family = "sorensen")

## Resultados
resultado_PA$beta.sor

## Data frame com os resultados
data.frame_PA <- data.frame(round(as.numeric(resultado_PA$beta.sor), 2),
                            round(as.numeric(resultado_PA$beta.sim), 2),
                            round(as.numeric(resultado_PA$beta.sne), 2))
colnames(data.frame_PA) <- c("Sorensen", "Simpson", "Aninhamento")
head(data.frame_PA)

## Dissimilaridade
prec_dis <- vegdist(precipitacao, method = "euclidian")
dados_prec <- as.numeric(prec_dis) 

## Juntando num Data frame
dados_dis <- data.frame(dados_prec, data.frame_PA)
head(dados_dis)

## ANOVA
# Avaliar a relação entre os valores de diversidade beta total (Sørensen) e precipitação.
anova_sore <-lm(Sorensen ~ dados_prec, data = dados_dis)
anova(anova_sore)


# Avaliar a relação entre os valores do componente substituição (Simpson) e precipitação
anova_simp <-lm(Simpson ~ dados_prec, data = dados_dis)
anova(anova_simp)

# Avaliar a relação entre os valores do componente aninhamento e precipitação
anova_anin <-lm(Aninhamento ~ dados_prec, data = dados_dis)
anova(anova_anin)

## Gráfico
ggplot(data = dados_dis, aes(x = dados_prec, y = Aninhamento)) + 
  geom_point(size = 4, shape = 21, fill = "darkorange") +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Diferença precipitação (mm)", 
       y = "Componente aninhamento da\n diversidade beta") +
  tema_livro()

## Interpretação dos resultados
## Diversidade beta para abundância
resultado_AB <- beta.pair.abund(composicao_especies, index.family = "bray")

## Data frame
# Vamos montar um data.frame com os resultados.
data.frame_AB <- data.frame(round(as.numeric(resultado_AB$beta.bray), 2),
                            round(as.numeric(resultado_AB$beta.bray.bal), 2),
                            round(as.numeric(resultado_AB$beta.bray.gra), 2))
colnames(data.frame_AB) <- c("Bray", "Balanceada", "Gradiente")
head(data.frame_AB)

## Agora vamos juntar os resultados com a precipitação
dados_dis_AB <- data.frame(dados_prec, data.frame_AB)

## ANOVA
# Avaliar a relação entre os valores de diversidade beta total e precipitação
anova_dis_AB <- lm(Bray ~ dados_prec, data = dados_dis_AB)
anova(anova_dis_AB)


# Avaliar a relação entre os valores do componente balanceada e precipitação
anova_balan <- lm(Balanceada ~ dados_prec, data = dados_dis_AB)
anova(anova_balan)
#> Analysis of Variance Table

# Avaliar a relação entre os valores do componente gradiente e precipitação
anova_grad <- lm(Gradiente ~ dados_prec, data = dados_dis_AB)
anova(anova_grad)

## gráfico para a variação balanceada da diversidade beta 
ggplot(data = dados_dis_AB, aes(x = dados_prec, y = Balanceada)) + 
  geom_point(size = 4, shape = 21, fill = "darkorange") +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Diferença precipitação (mm)",
       y = "Componente variação balanceada\n da diversidade beta") +
  tema_livro() 

## gráfico para a variação gradiente da diversidade beta
ggplot(data = dados_dis_AB, aes(x = dados_prec, y = Gradiente)) + 
  geom_point(size = 4, shape = 21, fill = "darkorange") +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Diferença precipitação anual (mm)", 
       y = "Componente gradiente de abundância\n da diversidade beta") +
  tema_livro()

