# Pacotes 
library(ade4)
library(ecodados)
library(tidyverse)
library(vegan) 
library(pvclust)
library(BiodiversityR)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(ape)
library(FactoMineR)
library(factoextra)
library(FD)
library(palmerpenguins)
library(GGally)
library(fields)
library(ade4)
library(ggord)
library(udunits2)
library(adespatial)
library(spdep)
library(mvabund)
library(reshape)

# Dados
sp_compos        <- ecodados::bocaina
species          <- ecodados::com_birds
env              <- ecodados::env_birds
xy               <- ecodados::birds.xy
bocaina.env      <- ecodados::bocaina.env
bocaina.xy       <- ecodados::bocaina.xy
anuros_permanova <- ecodados::anuros_permanova
macroinv         <- ecodados::macroinv
fish_comm        <- ecodados::fish_comm
data(mite)
data(doubs)
data(mite.env)

## Traduzir nomes para português 
colnames(penguins) <- c("especies", "ilha", "comprimento_bico", "profundidade_bico", "comprimento_nadadeira", "massa_corporal", "sexo", "ano")


# Análises de agrupamento
## Agrupamento hierárquico
### Matriz de similaridade com o coeficiente de Morisita-Horn
distBocaina <- vegdist(x = sp_compos, method = "horn")

## Agrupamento com a função hclust e o método UPGMA
dendro <- hclust(d = distBocaina, method = "average")

## Visualizar os resultados
plot(dendro, main = "Dendrograma", 
     ylab = "Similaridade (índice de Horn)",
     xlab="", sub="")

## Coeficiente de correlação cofenética
cofresult <- cophenetic(dendro)
cor(cofresult, distBocaina)

## Gráfico
plot(dendro, main = "Dendrograma", 
     ylab = "Similaridade (índice de Horn)",
     xlab="", sub="")
k <- 4
n <- ncol(sp_compos)
MidPoint <- (dendro$height[n-k] + dendro$height[n-k+1]) / 2
abline(h = MidPoint, lty=2)


### Outro exemplo

## Passo 1: transformar para distância de Chord
bocaina_transf <- disttransform(t(sp_compos), "chord")

## Passo 2: realizar pvclust com método average e distância euclidiana
analise <- pvclust(bocaina_transf, method.hclust = "average", method.dist = "euclidean", quiet = TRUE) 

## Passo 3: dendrograma
plot(analise, hang=-1, main = "Dendrograma com valores de P", 
     ylab = "Distância Euclideana",
     xlab="", sub="")
pvrect(analise)


### Análises

## Passo 1: transformar dados com Hellinger
bocaina_transf2 <- disttransform(t(bocaina), "hellinger")

## Passo 2: realizar pvclust com método average e distância euclidiana
analise2 <- pvclust(bocaina_transf2, method.hclust="average", method.dist="euclidean", quiet = TRUE) 

## Passo 3: dendrograma
plot(analise2, hang=-1, main = "Dendrograma com valores de P", 
     ylab = "Distância Euclideana",
     xlab="", sub="")
k <- 4
n <- ncol(sp_compos)
MidPoint <- (dendro$height[n-k] + dendro$height[n-k+1]) / 2
abline(h = MidPoint, lty=2)
pvrect(analise2)


## Agrupamento não-hierárquico (K-means)
## Verificar se existem localidades sem nenhuma ocorrência
rowSums(doubs$fish)

## Retirar a linha 8 (rio sem nenhuma ocorrência de peixe)
spe <- doubs$fish[-8,]

## Função do pacote vegan para normalizar os dados
spe.norm <- decostand(x = spe, method = "normalize") 

## K-Means
spe.kmeans <- kmeans(x = spe.norm, centers = 4, nstart = 100)
spe.kmeans

## Repetindo o K-Means
spe.KM.cascade <- cascadeKM(spe.norm, inf.gr = 2, sup.gr = 10, iter = 100, criterion = "ssi") 

## Resumo dos resultados
spe.KM.cascade$results
#>      2 groups  3 groups  4 groups  5 groups   6 groups  7 groups  8 groups  9 groups 10 groups
#> SSE 8.2149405 6.4768108 5.0719796 4.3015573 3.58561200 2.9523667 2.4840549 2.0521888 1.7599292
#> ssi 0.1312111 0.1685126 0.1409061 0.1299662 0.08693436 0.1481826 0.1267918 0.1134307 0.1226392

## Gráfico
plot(spe.KM.cascade, sortg = TRUE)


### Espécies indicadoras
fitofis <- c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4), rep(5, 4))

## Análise de espécies indicadoras
res_indval <- indval(t(sp_compos), fitofis)

# A função summary só exibe o resultado para as espécies indicadoras
summary(res_indval)

## Resultados
tab_indval <- cbind.data.frame(maxcls = res_indval$maxcls,
                               ind.value = res_indval$indcls,
                               P = res_indval$pval)
tab_indval

## Espécies
tab_indval[tab_indval$P < 0.05, ]

# Análises de Ordenação
## Ordenação irrestrita
### Análise de Componentes Principais (PCA)
#### Dados
aranhas <- data.frame(
  sp1 = c(5, 7, 2, 0, 0, 0, 0, 0),
  sp2 = c(0, 6, 3, 4, 0, 0, 0, 0),
  sp3 = c(0, 0, 0, 9, 12, 3, 0, 0),
  sp4 = c(0, 0, 0, 0, 4, 10, 8, 0),
  sp5 = c(0, 0, 0, 0, 0, 6, 9, 12),
  row.names = paste0("cidade", 1:8))

#### Centralização
aranha.cent <- as.data.frame(base::scale(aranhas, center = TRUE, scale=FALSE))

#### Matriz de covaiância
matriz_cov <- cov(aranha.cent)

#### Autovalores e autovetores
eigen_aranhas <- eigen(matriz_cov)
autovalores <- eigen_aranhas$values
autovetores <- as.data.frame(eigen_aranhas$vectors)
autovalores # eigenvalue

colnames(autovetores) <- paste("PC", 1:5, sep="")
rownames(autovetores) <- colnames(aranhas)
autovetores

#### Componentes principais
matriz_F <- as.data.frame(as.matrix(aranha.cent) %*% as.matrix(autovetores))
matriz_F

#### Porcentagem de explicação de cada eixo
100 * (autovalores/sum(autovalores))

#### Gráfico
ggplot(matriz_F, aes(x = PC1, y = PC2, label = rownames(matriz_F))) +
  geom_label() + 
  geom_hline(yintercept = 0, linetype=2) +
  geom_vline(xintercept = 0, linetype=2) +
  xlim(-10, 10) +
  tema_livro()

#### Verificar se existem NAs nos dados
sum(is.na(penguins))

#### Remover dados ausentes (NA), quando houver
penguins <- na.omit(penguins)

#### Manter somentes dados contínuos que pretende aplicar a PCA
penguins_trait <- penguins[, 3:6]

#### Compare com este código a variância das variáveis
penguins_trait %>% 
  dplyr::summarise(across(where(is.numeric), 
                          ~var(.x, na.rm = TRUE)))

#### Agora, veja o mesmo cálculo se fizer a padronização (scale.unit da função PCA)
penguins_pad <- decostand(x = penguins_trait, method = "standardize")
penguins_pad %>% 
  dplyr::summarise(across(where(is.numeric), 
                          ~var(.x, na.rm = TRUE)))
pca.p <- PCA(X = penguins_trait, scale.unit = TRUE, graph = FALSE)

#### Autovalores: porcentagem de explicação para usar no gráfico
pca.p$eig 

#### Visualização da porcentagem de explicação de cada eixo
# nota: é necessário ficar atento ao valor máximo do eixo 1 da análise para determinar o valor do ylim (neste caso, colocamos que o eixo varia de 0 a 70).
fviz_screeplot(pca.p, addlabels = TRUE, ylim = c(0, 70), main = "", 
               xlab = "Dimensões",
               ylab = "Porcentagem de variância explicada") 

#### Outros valores importantes
var_env <- get_pca_var(pca.p)

#### Escores (posição) das variáveis em cada eixo
var_env$coord 

#### Contribuição (%) das variáveis para cada eixo
var_env$contrib 

#### Loadings - correlação das variáveis com os eixos
var_env$cor 

#### Qualidade da representação da variável. Esse valor é obtido multiplicado var_env$coord por var_env$coord
var_env$cos2

#### Escores (posição) das localidades ("site scores") em cada eixo 
ind_env <- get_pca_ind(pca.p)

#### Variáveis mais importantes para o Eixo 1
dimdesc(pca.p)$Dim.1 

#### Variáveis mais importantes para o Eixo 2
dimdesc(pca.p)$Dim.2 

#### Gráfico de comparação morfológica
fviz_pca_biplot(X = pca.p, 
                geom.ind = "point", 
                fill.ind = penguins$especies, 
                col.ind = "black",
                alpha.ind = 0.7,
                pointshape = 21, 
                pointsize = 4,
                palette = c("darkorange", "darkorchid", "cyan4"),
                col.var = "black",
                invisible = "quali",
                title = NULL) +
  labs(x = "PC1 (68.63%)", y = "PC2 (19.45%)") + 
  xlim(c(-4, 5)) +
  ylim(c(-3, 3)) +
  tema_livro()

### Análises de Coordenadas Principais (PCoA)

#### Padronização dos dados com Hellinger
mite.hel <- decostand(x = mite, method = "hellinger") 

#### Cálculo da matriz de distância com método Bray-Curtis
sps.dis <- vegdist(x = mite.hel, method = "bray") 

#### PCoA
pcoa.sps <- pcoa(D = sps.dis, correction = "cailliez")

### Porcentagem de explicação do Eixo 1
100 * (pcoa.sps$values[, 1]/pcoa.sps$trace)[1]

### Porcentagem de explicação dos Eixo 2
100 * (pcoa.sps$values[, 1]/pcoa.sps$trace)[2]

### Porcentagem de explicação acumulada dos dois primeiros eixos 
sum(100 * (pcoa.sps$values[, 1]/pcoa.sps$trace)[1:2])

### Vizualização com GGPLOT

## Escores dos dois primeiros eixos
eixos <- pcoa.sps$vectors[, 1:2] 

## Combinar dados dos escores com um dado categórico de interesse para nossa pergunta
pcoa.dat <- data.frame(topografia = mite.env$Topo, eixos)

### Gráfico biplot da PCoA
ggplot(pcoa.dat, aes(x = Axis.1, y = Axis.2, fill = topografia, 
                     color = topografia, shape = topografia)) +
  geom_point(size = 4, alpha = 0.7) + 
  scale_shape_manual(values = c(21, 22)) + 
  scale_color_manual(values = c("black", "black")) + 
  scale_fill_manual(values = c("darkorange", "cyan4")) + 
  labs(x = "PCO 1 (49.11%)", y = "PCO 2 (14.30%)") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  tema_livro()

## Regressão de Componentes Principais (PCR)
### Dados
env_cont <- env[,-8]
env.pca <- PCA(env_cont, scale.unit = TRUE, graph = FALSE)
var_env <- get_pca_var(env.pca) 

### Contribuição (%) das variáveis para cada eixo
var_env$contrib 

### Loadings - correlação das variáveis com os eixos
var_env$cor 
ind_env <- get_pca_ind(env.pca)
env.pca$eig 

## Passo 1: obter os primeiros eixos 
pred.env <- ind_env$coord[, 1:3] 

## Passo 2: calcular a riqueza de espécies
riqueza <- specnumber(species)

## Passo 3: combinar os dois valores em um único data.frame
dat <- data.frame(pred.env, riqueza) 

## Regressão múltipla
mod1 <- lm(riqueza ~ Dim.1 + Dim.2 + Dim.3, data = dat)
par(mfrow = c(2, 2))
plot(mod1) # verificar pressupostos dos modelos lineares
summary(mod1) # resultados do  teste
dimdesc(env.pca)$Dim.1 

## Vizualização com ggplot
ggplot(dat, aes(x = Dim.1, y = riqueza)) + 
  geom_smooth(method = lm, fill = "#525252", color = "black") + 
  geom_point(size = 4, shape = 21, alpha = 0.7, color = "#1a1a1a", fill = "cyan4") +
  labs(x = "Gradiente ambiental (PC1)", y = "Riqueza de aves") + 
  tema_livro()


###  Outro exemplo
## Matriz de distância
env.dist <- gowdis(mite.env)

## PCoA
env.mite.pco <- pcoa(env.dist, correction = "cailliez")

## Porcentagem de explicação do Eixo 1
100 * (env.mite.pco$values[, 1]/env.mite.pco$trace)[1]

## Porcentagem de explicação dos Eixo 2
100 * (env.mite.pco$values[, 1]/env.mite.pco$trace)[2]

## Selecionar os dois primeiros eixos
pred.scores.mite <- env.mite.pco$vectors[, 1:2] 

## Juntar com os dados da área para fazer a figura 
mite.riqueza <- specnumber(mite)
pred.vars <- data.frame(riqueza = mite.riqueza, pred.scores.mite)

### Regressão múltipla 
mod.mite <- lm(riqueza ~ Axis.1 + Axis.2, data = pred.vars)
par(mfrow = c(2, 2))
plot(mod.mite)
summary(mod.mite)

### Vizualização com ggplot
g_acari_axi1 <- ggplot(pred.vars, aes(x = Axis.1, y = riqueza)) + 
  geom_smooth(method = lm, fill = "#525252", color = "black") + 
  geom_point(size = 4, shape = 21, alpha = 0.7, color = "#1a1a1a", fill="cyan4") + 
  labs(x = "Gradiente ambiental (PC1)", y = "Riqueza de ácaros") + 
  tema_livro()

g_acari_axi2 <- ggplot(pred.vars, aes(x = Axis.2, y = riqueza)) + 
  geom_smooth(method = lm, fill = "#525252", color = "black") + 
  geom_point(size = 4, shape = 21, alpha = 0.7, color = "#1a1a1a", fill = "darkorange") + 
  labs(x = "Gradiente ambiental (PC2)", y = "Riqueza de ácaros") + 
  tema_livro()

## Função para combinar os dois plots em uma única janela
grid.arrange(g_acari_axi1, g_acari_axi2, nrow = 1)


# Ordenação restrita
## Análise de Redundância (RDA) - Exemplo 1
### Passo 1: transformação de hellinger da matriz de espécies
### caso tenha dados de abundância.
species.hel <- decostand(x = species, method = "hellinger")

## Passo 2: selecionar variáveis importantes
# Para isso, é necessário remover a variável categórica.
env.contin <- env[, -8]

## Evite usar variáveis muito correlacionadas
sel.vars <- forward.sel(species.hel, env.contin)
sel.vars$variables
env.sel <- env[,sel.vars$variables]

## Passo 3: padronizar matriz ambiental (somente variáveis contínuas)
env.pad <- decostand(x = env.sel, method = "standardize")

## Matriz final com variáveis preditoras
env.pad.cat <- data.frame(env.pad, altitude = env$altitude)

## RDA com dados selecionados e padronizados
rda.bird <- rda(species.hel ~ rain.jul + maxi.jul + altitude, data = env.pad.cat)

## Para interpretar, é necessário saber a significância dos eixos para representar a relação entre as variáveis preditoras e a composição de espécies
res.axis <- anova.cca(rda.bird, by = "axis") 
res.axis

## Em seguida, é possível identificar quais são as variáveis que contribuem ou que mais contribuem para a variação na composição de espécies 
res.var <- anova.cca(rda.bird, by = "term") ## Qual variável?
res.var

## Além disso, é possível obter o valor do R2 do modelo
r_quadr <- RsquareAdj(rda.bird)
r_quadr

## Ordenação multi-escala (MSO) para entender os resultados da ordenação em relação à distância geográfica
bird.rda <- mso(rda.bird, xy, grain = 1, permutations = 99)
msoplot(bird.rda)

## Triplot da RDA
ggord(rda.bird, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + 
  tema_livro()

## Análise de Redundância parcial (RDAp)
## Passo 1: Gerar um arquivo LIST W: list binária de vizinhança
mat_knn <- knearneigh(as.matrix(xy), k = 2, longlat = FALSE)
mat_nb <- knn2nb(mat_knn, sym = TRUE)
mat_listw <- nb2listw(mat_nb, style = "W")
mat_listw

## Passo 2: Listar os métodos "candidatos" para obter a matriz SWM
MEM_mat <- scores.listw(mat_listw, MEM.autocor = "positive")
candidates <- listw.candidates(xy, nb = c("gab", "mst", "dnear"), 
                               weights = c("binary", "flin"))

## Passo 3: Selecionar a melhor matriz SWM e executar o MEM
W_sel_mat <- listw.select(species.hel, candidates, MEM.autocor = "positive",
                          p.adjust = TRUE, method = "FWD")

## Passo 4: Matriz dos preditores espaciais escolhidos (MEMs)
spatial.pred <- as.data.frame(W_sel_mat$best$MEM.select)

## É necessário atribuir os nomes das linhas
rownames(spatial.pred) <- rownames(xy) 

## Combinar variáveis ambientais e espaciais em um único data.frame
pred.vars <- data.frame(env.pad.cat, spatial.pred)

## RDA parcial
rda.p <- rda(species.hel ~
               rain.jul + maxi.jul + altitude + # Preditores ambientais
               Condition(MEM1 + MEM2 + MEM4 + MEM5), # Preditores espaciais
             data = pred.vars)

## Interpretação
# Para interpretar, é necessário saber a significância dos eixos para representar a relação entre as variáveis preditoras e a composição de espécies.
res.p.axis <- anova.cca(rda.p, by = "axis") 
res.p.axis

## Contribuição
# Em seguida, é possível identificar quais são as variáveis que contribuem ou que mais contribuem para a variação na composição de espécies.
res.p.var <- anova.cca(rda.p, by = "term") ## Qual variável?
res.p.var
RsquareAdj(rda.p)

## Padrão espacial na composição de espécies
pca.comp <- dudi.pca(species.hel, scale = FALSE, scannf = FALSE)
moran.comp <- moran.mc(pca.comp$li[, 1], mat_listw, 999)

## Padrão espacial das variáveis ambientais
env$altitude <- as.factor(env$altitude)
ca.env <- dudi.hillsmith(env, scannf = FALSE)
moran.env <- moran.mc(ca.env$li[, 1], mat_listw, 999)

## Estrutura espacial na composição de espécies?
moran.comp

## Estrutura espacial na variação ambiental?
moran.env

## Triplot da RDA
# Partição de variância.
pv.birds <- varpart(species.hel, env.pad.cat, spatial.pred)
plot(pv.birds, cutoff = -Inf)
mtext("Diagrama de Venn: partição de variância")

## Análise de Redundância baseada em distância (db-RDA)
## Passo 1: transformação de hellinger da matriz de espécies
# caso tenha dados de abundância.
species.hel <- decostand(species, "hellinger")

## Passo 2: gerar matriz de distância (Diversidade beta: Bray-Curtis)
dbeta <- vegdist(species.hel, "bray")

## Passo 2: dbRDA 
dbrda.bird <- capscale (dbeta~rain.jul+maxi.jul+altitude, data=env.pad.cat)

# Para interpretar, é necessário saber a significância dos eixos para representar a relação entre as variáveis preditoras e a composição de espécies
db.res.axis <- anova.cca(dbrda.bird, by = "axis") 
db.res.axis

# Em seguida, é possível identificar quais são as variáveis que contribuem ou que mais contribuem para a variação na composição de espécies 
db.res.var <- anova.cca(dbrda.bird, by = "term") ## Qual variável?
db.res.var

# Além disso, é possível obter o valor do R2 do modelo
db_r_quadr <- RsquareAdj(dbrda.bird)
db_r_quadr


# PERMANOVA
## Composição de espécies padronizar com método de Hellinger
species.hel <- decostand(x = species, method = "hellinger")

## Matriz de distância com método Bray-Curtis
sps.dis <- vegdist(x = species.hel, method = "bray")

## Verifica correlação entre as variáveis
ggpairs(env) +
  tema_livro()

## PERMANOVA
perm.aves <- adonis2(sps.dis ~ mini.jan + rain.tot + altitude, data = env2)
perm.aves ### Diferenças entre os tratamentos?

## BETADISPER
betad.aves <- betadisper(sps.dis, env2$altitude)
permutest(betad.aves)


## Após verificar a estrutura de correlação, vamos manter somente três variáveis
env2 <- env[, c("mini.jan", "rain.tot", "altitude")]

## Matriz de distância representando a variação na composição de espécies (método Bray-Curtis)
as.matrix(sps.dis)[1:6, 1:6]

## É preciso calcular uma primeira "melhor" solução do nMDS
sol1 <- metaMDS(sps.dis)

## Melhor solução
# Depois, executar a mesma função, mas utilizando uma "melhor solução inicial" para evitar resultdos subótimos no nMDS .
nmds.beta <- metaMDS(sps.dis, previous.best = sol1)

## Stress
# O stress é o valor mais importante para interpretar a qualidade da ordenação.
nmds.beta$stress # valor ideal entre 0 e 0.2

## Exportar os valores para fazer gráfico 
dat.graf <- data.frame(nmds.beta$points, altitude = env2$altitude)

## Definir os grupos ("HULL") para serem categorizados no gráfico 
grp.mon <- dat.graf[dat.graf$altitude == "Montanhoso", ][chull(dat.graf[dat.graf$altitude == "Montanhoso", c("MDS1", "MDS2")]), ]

grp.int <- dat.graf[dat.graf$altitude == "Intermediário", ][chull(dat.graf[dat.graf$altitude == "Intermediário", c("MDS1", "MDS2")]), ]

grp.pla <- dat.graf[dat.graf$altitude == "Plano", ][chull(dat.graf[dat.graf$altitude == "Plano", c("MDS1", "MDS2")]), ]

## Combinar dados dos grupos para cada Convex Hull
hull.data <- rbind(grp.mon, grp.int, grp.pla) 

## Gráfico
ggplot(dat.graf, aes(x = MDS1, y = MDS2, color = altitude, shape = altitude)) + 
  geom_point(size = 4, alpha = 0.7) + 
  geom_polygon(data = hull.data, aes(fill = altitude, group = altitude), alpha = 0.3) + 
  scale_color_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  scale_fill_manual(values = c("darkorange", "darkorchid", "cyan4")) +
  labs(x = "NMDS1", y = "NMDS2") + 
  tema_livro()

# Mantel e Mantel Parcial
## Matriz de distância geográfica
Dist.km <- as.dist(rdist.earth(bocaina.xy, miles=F)) # matriz de distância geográfica (geodésica) considerando a curvatura da Terra

## Padronizações
comp.bo.pad <- decostand(t(bocaina), "hellinger")
env.bocaina <- decostand(bocaina.env[,-9], "standardize")

## Dissimilaridades
dissimil.bocaina <- vegdist(comp.bo.pad, "bray")
dissimil.env <- vegdist(env.bocaina, "euclidian")

## Mantel
# Espaço vs. ambiente
mantel(Dist.km, dissimil.env) 

## Preparar os dados para o gráfico
matrix.dist.env <- data.frame(x = melt(as.matrix(Dist.km))$value, 
                              y = melt(as.matrix(dissimil.env))$value)
## Gráfico
ggplot(matrix.dist.env , aes(x, y)) +
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  labs(x = "Distância geográfica (km)", 
       y = "Dissimilaridade Ambiental") +
  tema_livro()                              

## Mantel
# Espaço vs. composição
mantel(Dist.km, dissimil.bocaina) 


## Preparar os dados para o gráfico
matrix.dist.bocaina <- data.frame(x = melt(as.matrix(Dist.km))$value, 
                                  y = melt(as.matrix(dissimil.bocaina))$value)

## Gráfico
ggplot(matrix.dist.bocaina , aes(x, y)) +
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  labs(x = "Distância geográfica (km)", 
       y = "Dissimilaridade (Bray-Curtis)") +
  tema_livro()

## Mantel
# Ambiente vs. composição
mantel(dissimil.env, dissimil.bocaina) 

## Preparar os dados para o gráfico
matrix.bocaina.env <- data.frame(x = melt(as.matrix(dissimil.bocaina))$value, 
                                 y = melt(as.matrix(dissimil.env))$value)

## Gráfico
ggplot(matrix.bocaina.env, aes(x, y)) +
  geom_point(size = 4, shape = 21, fill = "darkorange", alpha = 0.7) +
  labs(x = "Similaridade Bocaina", 
       y = "Similaridade Ambiental") +
  tema_livro()

### Interpretação dos resultados
## Mantel Parcial
mantel.partial(dissimil.bocaina, dissimil.env, Dist.km) #comp. vs. ambiente controlando o espaço

# Mantel espacial com modelo nulo restrito considerando autocorrelação espacial
## Mantel
compos_espac <- mantel.randtest(sqrt(dissimil.bocaina), Dist.km)
compos_espac

## Minimum Spanning Tree
nb.boc <- mst.nb(Dist.km) # calcula uma Minimum Spanning Tree
plot(nb.boc, bocaina.xy)

## Conversão
lw <- nb2listw(nb.boc)

## Mantel espacial
msr(compos_espac, lw, method = "pair")

## Mantel comum
compos_espac


# PROCRUSTES e PROTEST
## Fixar a amostragem
set.seed(1001) 

## Matrizes de distância de PCoA 
d_macro <- vegdist(macroinv, "bray")
pcoa_macro <- cmdscale(d_macro)

d_fish <- vegdist(fish_comm, "bray")
pcoa_fish <- cmdscale(d_fish)

## PROCRUSTES
concord <- procrustes(pcoa_macro, pcoa_fish)
protest(pcoa_macro, pcoa_fish)

## Gráfico
plot(concord, main = "", 
     xlab = "Dimensão 1",
     ylab = "Dimensão 2")

# Métodos multivariados baseados em modelos
#retirando a coluna que contém o fator e deixando apenas dados de abundância
anuros_abund <- as.data.frame(anuros_permanova[,-28])

#incluindo apenas o fator a ser testado no modelo
grupos <- factor(anuros_permanova[,28])

## Média-variância
# criando um objeto da classe mvabund com os dados de abundância
abund_tr <- mvabund(anuros_abund) 
meanvar.plot(abund_tr)

# Exemplo 1
## Gráfico
plot(abund_tr ~ grupos, cex.axis = 0.8, cex = 0.8)
#> 
#>  PIPING TO 2nd MVFACTOR

## Modelo
modelo1 <- manyglm(abund_tr ~ grupos, family = "negative.binomial")

## Diagnose
plot(modelo1)

## Resultados
summary(modelo1)
### ANOVA
anova(modelo1, p.uni = "adjusted")


