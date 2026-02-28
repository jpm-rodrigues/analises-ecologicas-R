## Pacotes
library(FD)
library(ade4)
library(ecodados)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(picante)
library(vegan)
library(SYNCSA)
library(GGally)
library(FD)
library(betapart)
library(nlme)
library(ape)
library(TPD)
library(cati)
library(kableExtra)

## Dados e funções
comun_fren_dat <- ecodados::fundiv_frenette2012a_comu
ambie_fren_dat <- ecodados::fundiv_frenette2012a_amb
trait_fren_dat <- ecodados::fundiv_frenette2012a_trait
trait_dat <- ecodados::fundiv_barbaro2009a_trait
comun_dat <- ecodados::fundiv_barbaro2009a_comu
ambie_dat <- ecodados::fundiv_barbaro2009a_amb
trait_baselga <- ecodados::trait_baselga
comm_baselga <- ecodados::comm_baselga
anuros_comm <- ecodados::anuros_comm
traits <- ecodados::traits
env <- ecodados::env

# Definindo a dis(similaridade) entre espécies
## Exemplo1: variáveis contínuas
## PCoA dos atributos contínuos

# 1. Padronização dos dados
trait_pad <- decostand(trait_fren_dat, "standardize")
euclid_dis <- vegdist(trait_pad, "euclidean")

# 2. PCoA
# Resultados são idênticos aos resultados de uma PCA.
pcoa_traits_cont <- pcoa(euclid_dis, correction = "cailliez")

# 3. Exportandos dados para gráfico
# Ao usar '[,1:2]' você irá selecionar os dois primeiros eixos.
eixos_cont <- as.data.frame(pcoa_traits_cont$vectors[, 1:2])

# 4. Gráfico de ordenação
plot_trait_cont <- ggplot(eixos_cont, aes(x = Axis.1, y = Axis.2)) +
  geom_point(pch = 21, size = 4, color = "black", alpha = 0.7, fill = "red2") +
  geom_text_repel(aes(Axis.1, Axis.2, label = rownames(eixos_cont))) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "PCO 1", y = "PCO 2", title = "Dados contínuos") +
  tema_livro()
plot_trait_cont

## Exemplo 2: variáveis categóricas
## PCoA dos atributos categóricos

# 1. Selecionar somente os atributos categóricos
trait_cat <- trait_dat %>%
  dplyr::select_if(is.character)

# 2. Calcular a distância de Gower
dist_categ <- gowdis(trait_cat)

# 3. PCoA da matriz de distância funcional (Gower)
pcoa_traits_cat <- pcoa(dist_categ, correction = "cailliez")

# 4. Exportar dados (escores) para ggplot
eixos_cat <- as.data.frame(pcoa_traits_cat$vectors[, 1:2]) # Selecionar os dois primeiros eixos

# 5. Gráfico de ordenação
plot_trait_cat <- ggplot(eixos_cat, aes(x = Axis.1, y = Axis.2)) +
  geom_point(pch = 21, size = 4, alpha = 0.7, color = "black", fill = "cyan4") +
  geom_text_repel(aes(Axis.1, Axis.2, label = rownames(eixos_cat))) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "PCO 1", y = "PCO 2", title = "Dados categóricos") +
  tema_livro()
plot_trait_cat

## Exemplo 3: variáveis mistas
## PCoA dos atributos mistos

## 1. Verifique a classe de todos os traits e veja se estão de acordo com sua expectativa
trait_dat %>%
  dplyr::summarise_all(class) %>%
  tidyr::gather(variable, class)

## 2. Neste exemplo, algumas variáveis que são ordinais (regio e body)
# foram reconhecidas como numéricas ou categóricas.
trait_dat$regio <- as.ordered(trait_dat$regio)
trait_dat$body <- as.ordered(trait_dat$body)

## 3. Combinar cada conjunto de atributos de acordo com sua natureza em um
# data.frame separado.
# 3.1. Categóricos.
trait_categ <- cbind.data.frame(
  trend = trait_dat$trend,
  redlist = trait_dat$redlist,
  biog = trait_dat$biog,
  activ = trait_dat$activ,
  diet = trait_dat$diet,
  winter = trait_dat$winter,
  color = trait_dat$color,
  breed = trait_dat$breed,
  wing = trait_dat$wing,
  period = trait_dat$period
)

# 3.2 Ordinais.
trait_ord <- cbind.data.frame(
  regio = trait_dat$regio,
  body = trait_dat$body
)
rownames(trait_categ) <- rownames(trait_dat)
rownames(trait_ord) <- rownames(trait_dat)

# Agora, combinar os dois data.frames em uma lista chamada "ktab".
ktab_list <- ktab.list.df(list(trait_categ, trait_ord))

# Por fim, calcular a distância funcional entre as espécies.
# Em "type", a letra "N" indica variável categórica (ou nominal),
# enquanto a letra "O" indica variável ordinal.
dist_mist <- dist.ktab(ktab_list, type = c("N", "O"))

## Visualize os dados com uma PCoA
pcoa_traits_mist <- pcoa(dist_mist, correction = "cailliez")
eixos_mist <- as.data.frame(pcoa_traits_mist$vectors[, 1:2])

plot_trait_mist <- ggplot(eixos_mist, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    pch = 21, size = 4, alpha = 0.7,
    color = "black", fill = "darkorange"
  ) +
  geom_text_repel(aes(Axis.1, Axis.2, label = rownames(eixos_mist))) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "PCO 1", y = "PCO 2", title = "Dados mistos") +
  tema_livro()
plot_trait_mist

## combinar os dois gráficos (baseado em variáveis categóricas e em
## variáveis mistas) para comparar as duas medidas de distância

grid.arrange(plot_trait_cat, plot_trait_mist, ncol = 2)

# Métricas de diversidade funcional (alfa)
# Riqueza funcional
## Exemplo prático
## Estrutura dos dados
# matriz de distância: distância entre as seis primeiras espécies
as.matrix(dist_mist)[1:6, 1:6]

# composição de espécies: seis primeiras espécies nas seis primeiras localidades
head(comun_dat)[1:6, 1:6]

## Calculando a riqueza de espécies por comunidade
richness <- dbFD(dist_mist, comun_dat)$nbsp
head(richness)

## Functional Richness
fric <- dbFD(dist_mist, comun_dat)$FRic
head(fric)

### Functional Diversity
### Passo 1: análise de agrupamento para criar o dendrograma.
dend <- hclust(dist_mist, "average")

### Passo 2: transformar o dengrograma em um arquivo da classe phylo.
tree_dend <- as.phylo(dend)

### Passo 3: calcular o valor da diversidade funcional.
FD <- pd(comun_dat, tree_dend)$PD
head(FD)

# Divergência funcional
## "Functional Divergence" s
fdiv <- dbFD(dist_mist, comun_dat)$FDiv
head(fdiv)
## "Functional Dispersion"
fdis <- dbFD(dist_mist, comun_dat)$FDis
head(fdis)

# Regularidade funcional
## "Functional evenness"
feve <- dbFD(dist_mist, comun_dat)$FEve
head(feve)

# Correlação entre as métricas de diversidade funcional (alfa)
## Gráficos
## Você pode criar uma tabela com os resultados de todas as métricas
metricas <- data.frame(
  richness = richness,
  FD_gp = FD,
  fric = fric,
  fdiv = fdiv,
  fdis = fdis,
  feve = feve
)
head(metricas)
## Gráfico para comparar o comportamento das métricas
ggpairs(metricas) + tema_livro()

# Métricas de diversidade funcional (beta)
## Exemplo 4
## Partição da Diversidade beta (Método Baselga)
fun_beta_multi <- functional.beta.multi(
  x = comm_baselga,
  trait = trait_baselga, index = "jaccard"
)
fun_beta_multi
## Partição da Diversidade beta (Método Baselga)
fun_beta <- functional.beta.pair(
  x = comm_baselga,
  trait = trait_baselga, index = "jaccard"
)

# Os códigos abaixo permitem extrair a matriz de distância (par a par) com a partição em substituição e nestedness
fun_turnover <- fun_beta$funct.beta.jne
fun_nestedness <- fun_beta$funct.beta.jtu
fun_jaccard <- fun_beta$funct.beta.jac

## Gráfico de comparação do substituição e aninhamento
dat_betapart <- data.frame(
  turnover = as.numeric(fun_turnover),
  nested = as.numeric(fun_nestedness)
)

## Gráfico
plot_betapart <- ggplot(dat_betapart, aes(x = turnover, y = nested)) +
  geom_point(pch = 21, size = 4, alpha = 0.7, color = "black", fill = "#525252") +
  labs(x = "Beta Diveristy (Substituição)", y = "Beta Diveristy (Aninhamento)") +
  tema_livro()
plot_betapart

# Composição Funcional (Community Wegihed Means - CWM)
## Matriz T
head(trait_baselga)
## Matriz X
head(comm_baselga)
## Função functcomp calcula o cwm para combinar as matrizes T e X
cwm_ex <- functcomp(trait_baselga, as.matrix(comm_baselga))
cwm_ex

## Exemplo 5

## Passo 1: calcular a distância funcional
trait_pad <- decostand(trait_fren_dat, "standardize")
euclid_dis <- vegdist(trait_pad, "euclidean")

## Passo 2: calcular a Divergência funcional (FDis) e Regularidade Funcional (FEve)
fdis <- dbFD(euclid_dis, comun_fren_dat)$FDis # Fdis=0 em locais com somente uma espécie
feve <- dbFD(euclid_dis, comun_fren_dat)$FEve

## Passo 3: Utilizar um modelo linear para comparar o efeito da aridez sobre FDis (predição 1) e FEve (predição 2)
# Combinar dados em um data.frame.
lm_dat <- data.frame(aridez = ambie_fren_dat$Aridity, fdis = fdis, feve = feve)

# Modelo 1
mod1 <- lm(fdis ~ aridez, data = lm_dat)

# Conclusão: a aridez não tem efeito sobre a divergência funcional
anova(mod1)

# Modelo 2
mod2 <- lm(feve ~ aridez, data = lm_dat)

# Conclusão: a aridez não tem efeito sobre a regularidade funcional
anova(mod2)

## Passo 4: gráfico para visualizar os dois resultados
# Gráfico modelo 1.
plot_pred1 <- ggplot(lm_dat, aes(x = aridez, y = fdis)) +
  geom_point(pch = 21, size = 4, alpha = 0.7, color = "black", fill = "darkorange") +
  labs(x = "Aridez", y = "Divergência Funcional (FDis)") +
  tema_livro()

# Gráfico modelo 2.
plot_pred2 <- ggplot(lm_dat, aes(x = aridez, y = feve)) +
  geom_point(pch = 21, size = 4, alpha = 0.7, color = "black", fill = "cyan4") +
  labs(x = "Aridez", y = "Regularidade Funcional (FEve)") +
  tema_livro()

## Visualização dos dois gráficos em um única janela
grid.arrange(plot_pred1, plot_pred2, ncol = 2)


## Exemplo 6
## Passo 1: CWM
cwm_fren <- functcomp(trait_pad, as.matrix(comun_fren_dat))
head(cwm_fren)

## Passo 2: calcular a distância funcional
cwm_dis <- vegdist(cwm_fren, "euclidean")

## Passo 3: testar se a composição funcional varia entre as áreas com uma PERMANOVA
perman_fren <- adonis2(cwm_fren ~ Grazing, data = ambie_fren_dat)

## Passo 4: comparar a variação dentro de cada grupo com Betadisper
betad_fren <- betadisper(cwm_dis, ambie_fren_dat$Grazing)
permutest(betad_fren)

## Passo 5: PCoA

cwm_pcoa <- pcoa(D = cwm_dis, correction = "cailliez")
pcoa_eixos <- cwm_pcoa$vectors[, 1:2]
pcoa_dat <- data.frame(pastagem = ambie_fren_dat$Grazing, pcoa_eixos)

## Passo 6: definir os grupos ("HULL") para serem categorizados no gráfico

grp.Grazed <- pcoa_dat[pcoa_dat$pastagem == "Grazed", ][chull(pcoa_dat[pcoa_dat$pastagem == "Grazed", c("Axis.1", "Axis.2")]), ]
grp.Ungrazed <- pcoa_dat[pcoa_dat$pastagem == "Ungrazed", ][chull(pcoa_dat[pcoa_dat$pastagem == "Ungrazed", c("Axis.1", "Axis.2")]), ]
hull_cwm <- rbind(grp.Grazed, grp.Ungrazed)

## Passo 7: Gráfico biplot

100 * (cwm_pcoa$values[, 1] / cwm_pcoa$trace)[1] # % de explicação do eixo 1
100 * (cwm_pcoa$values[, 1] / cwm_pcoa$trace)[2] # % de explicação do eixo 2

ggplot(pcoa_dat, aes(x = Axis.1, y = Axis.2, color = pastagem, shape = pastagem)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_polygon(data = hull_cwm, aes(fill = pastagem, group = pastagem), alpha = 0.3) +
  scale_color_manual(values = c("darkorange", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "cyan4")) +
  labs(x = "PCO 1 (53.6%)", y = "PCO 2 (24.6%)") +
  tema_livro()

# Variação Intraespecífica
## Exemplo 7
## Pergunta 1
## Dados necessários
# Matriz de traits.
head(traits)

## Partição da variação intra e interespecífica
## Passo 1: Tamanho corporal
mod_body_size <- aov(body_size ~ Species, data = traits)
summary(mod_body_size)

# Contribuição da variação intra-específica para o tamanho corporal.
itv_BS <- 100 * (73.35 / (95.92 + 73.35))
itv_BS

## Passo 2: Biomassa
mod_biomass <- aov(biomass ~ Species, data = traits)
summary(mod_biomass)

# Contribuição da variação intra-específica para a biomassa.
itv_biomass <- 100 * (96.22 / (118.17 + 96.22))
itv_biomass

## Passo 3: Tamanho do olho
mod_eye_size <- aov(eye_size ~ Species, data = traits)
summary(mod_eye_size)

## Contribuição da variação intra-específica para o tamanho do olho.
itv_eye_size <- 100 * (78.39 / (203.09 + 78.39))
itv_eye_size

## Passo 4: Achatamento dorso-ventral
mod_flatness <- aov(flatness ~ Species, data = traits)
summary(mod_flatness)

## Contribuição da variação intra-específica para o achatamento dorso-ventral.
itv_flatness <- 100 * (22.13 / (104.48 + 22.13))
itv_flatness

## Passo 5: Combinar os valores de cada trait em um vetor
valores <- c(itv_BS, itv_biomass, itv_eye_size, itv_flatness)

# Passo 6: Combinar valores e traits em um data.frame.
itv_results <- data.frame(
  trait = c("body_size", "biomass", "eye_size", "flatness"),
  itv_explic = valores
)

## Tabela com resultados da explicação atribuida para a variação intraespecífica
itv_results %>%
  mutate("explained intraspecific variance" = round(itv_explic, 2)) %>%
  dplyr::select(trait, "explained intraspecific variance")

## Pergunta 2
## Dados necessários
# Matriz de traits sem nomes de espécies ou localidades
trait_m <- traits[, c("body_size", "biomass", "eye_size", "leg_size", "flatness")]
head(trait_m)
#>   body_size biomass eye_size leg_size flatness
#> 1     2.405   2.291    3.104    0.450    0.794
#> 2     1.882   2.039    2.926    0.345    1.063
#> 3     0.699   0.342    0.782    0.104    3.055
#> 4     0.725   0.598    1.120    0.136    2.759
#> 5     0.448   0.385    0.844    0.107    3.557
#> 6     0.640   0.470    0.861    0.093    3.420
trait_decomp <- decompCTRE(
  traits = trait_m, sp = traits$Species,
  ind.plot = traits$pond, print = FALSE
)
barplot.decompCTRE(trait_decomp)

## Pergunta 3

## Dados necessários
# Matriz de traits.
head(traits)

# Matriz de comunidades e padronização para abundância relativa
head(anuros_comm)
anuros_comm_rel <- decostand(anuros_comm, "total")

# Variáveis ambientais.
head(env)

## Prearação da matriz para receber os resultados do `for`
wITVResults <- data.frame(ITV = matrix(ncol = 1, nrow = length(unique(traits$pond))))
rownames(wITVResults) <- unique(traits$pond)

for (i in 1:length(unique(traits$pond))) {
  commAux <- subset(traits, traits$pond == unique(traits$pond)[i])
  commAux$Species <- droplevels(factor(commAux$Species))
  spNames <- unique(commAux$Species)
  relAbund <- anuros_comm_rel[i, as.character(spNames)]
  traitsVector <- commAux$body_size
  spVector <- commAux$Species
  wITVResults[i, 1] <- wITV(spIDs = spVector, traitVals = traitsVector, relAbund = relAbund)
}

wITVResults$ITV

## Remover NAs para executar o modelo linear
env2 <- na.omit(env)
head(env2)

## Modelo linear
mod_itv <- lm(wITVResults$ITV ~ depth + area + dits_bt_pond + dist_for, data = env)

## Testar pressuposto da análise
par(mfrow = c(2, 2))
plot(mod_itv)

## Resultado
summary(mod_itv)

