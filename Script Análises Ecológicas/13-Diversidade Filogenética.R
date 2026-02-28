## Pacotes
library(devtools)
library(ecodados)
library(V.PhyloMaker)
library(vegan)
library(ggplot2)
library(GGally)
library(ggpubr)
library(picante)
library(phytools)
library(ape)
library(geiger)
library(phyloregion)
library(pez)
library(reshape2)
library(betapart)

## Dados
minha_arvore <- ecodados::filogenia_aves
especies_plantas <- ecodados::sp_list
comunidade <- ecodados::comm
composicao_especies <- ecodados::composicao_aves_filogenetica
filogenia_aves <- ecodados::filogenia_aves
precipitacao <- precipitacao_filogenetica

# Manipulação de filogenias 
## Gráfico
plot.phylo(minha_arvore,
           type = "phylogram", show.tip.label = TRUE,
           show.node.label = TRUE, edge.color = "black", edge.width = 1.5, 
           tip.color = "black", cex = 0.45, label.offset = 2) 

## Gráfico
plot.phylo (minha_arvore, type = "fan", show.tip.label = TRUE, 
            show.node.label = TRUE, edge.color = "blue", edge.width = 1.5, 
            tip.color = "black", cex = 0.45, label.offset = 2) 

## Acessar informações da filogenia
## Nomes
names(minha_arvore)
### Nome das espécies
minha_arvore$tip.label

### Comprimento dos ramos
minha_arvore$edge.length

## Remover espécies da filogenia
# Vamos criar um novo nome para o objeto e excluir as espécies Leucopternis polionotus e Aramides saracura da filogenia
filogenia_cortada <- drop.tip(minha_arvore, c("Leucopternis_polionotus", "Aramides_saracura"))
filogenia_cortada

## Adicionar espécies à filogenia
# Vamos inserir as espécies Megascops_sp1, Carponis_sp, Strix_sp1, Strix_sp2 e
# Strix_sp3 na filogenia
Megascops <- c("Megascops_sp1")
Carpornis <- c("Carpornis_sp1")
Strix <- c("Strix_sp1", "Strix_sp2", "Strix_sp3")

# Inserindo espécies como politomias
filogenia_nova <- add.species.to.genus(force.ultrametric(minha_arvore, message = FALSE), Megascops)
filogenia_nova <- add.species.to.genus(force.ultrametric(filogenia_nova, message = FALSE), Carpornis)

## Adicionar várias espécies à filogenia
# Para inserir mais de uma espécie dentro do gênero, vamos utilizar um loop.  
for(i in 1:length(Strix)) 
  filogenia_nova <- add.species.to.genus(force.ultrametric(filogenia_nova, message = FALSE),
                                         Strix[i], where = "root")
## Gráfico
plot(filogenia_nova, cex = 0.5, no.margin = TRUE)

## phylo.maker #############
# A função phylo.maker usa uma filogenia default de plantas (i.e. GBOTB.extended).
# Caso você queira utilizar outra filogenia, é só alterar o argumento tree
novas_filogenias <- phylo.maker(especies_plantas,
                                tree = GBOTB.extended,
                                scenarios = c("S1","S2","S3"))
## Gráfico
par(mfrow = c(1, 2))
plot.phylo(novas_filogenias$scenario.1, cex = 0.5, main = "Cenário 1")
plot.phylo(novas_filogenias$scenario.3, cex = 0.5, main = "Cenário 3")
dev.off()

# Métricas de diversidade alfa filogenética
# Riqueza da diversidade alfa filogenética
# Phylogenetic diversity (PD)

## Conferir os nomes das espécies
name.check(filogenia_aves, t(composicao_especies))

## Colocar os nomes das espécies do data frame na mesma ordem que aparecem na filogenia
composicao_especies_P <- match.phylo.comm(phy = filogenia_aves, comm = composicao_especies)$comm

## Phylogenetic diversity (PD)
## Calculando a métrica de diversidade filogenética proposta por Faith (1992).
resultados_PD <- pd(composicao_especies_P, filogenia_aves)

## Mostra o valor de PD e riqueza de espécies para cada comunidade.
resultados_PD

# Phylogenetic Species Richness (PSR)
resultados_PSR <- psr(composicao_especies_P,filogenia_aves)
resultados_PSR 

# Phylogenetic Endemism (PE)
## Transformando data.frame em matriz.
dados_matriz <- as.matrix(composicao_especies_P)

## Análise.
resultados_PE <- phylo_endemism(dados_matriz, filogenia_aves, 
                                weighted = TRUE)

## Mostra os valores de PE para cada comunidade.
resultados_PE 

# Species Evolutionary Distinctiveness (ED)
## Análise.
resultados_ED <- evol.distinct(filogenia_aves)

## Mostra os valores de ED para cada espécie.
head(resultados_ED)

# Divergência da diversidade alfa filogenética
# Mean Pairwise Distance (MPD)

## Análise com dados de incidência das espécies nas comunidades.
resultados_MPD_PA <- mpd(composicao_especies_P, cophenetic(filogenia_aves), 
                         abundance.weighted = FALSE)

## Mostra os valores de MPD para cada comunidade.
resultados_MPD_PA 

## Análise com dados de abundância das espécies nas comunidades.
resultados_MPD_AB <- mpd(composicao_especies_P, cophenetic(filogenia_aves), 
                         abundance.weighted = TRUE)

## Mostra os valores de MPD para cada comunidade.
resultados_MPD_AB 

# Mean Nearest Taxon Distance (MNTD)
## Análise com dados de presença e ausência das espécies nas comunidades.
resultados_MNTD_PA <- mntd(composicao_especies_P, cophenetic(filogenia_aves), 
                           abundance.weighted = FALSE)

## Mostra os valores de MPD para cada comunidade.
resultados_MNTD_PA 

## Análise com dados de abundância das espécies nas comunidades.
resultados_MNTD_AB <- mntd(composicao_especies_P, cophenetic(filogenia_aves), 
                           abundance.weighted = TRUE)

## Mostra os valores de MPD para cada comunidade.
resultados_MNTD_AB 

# Phylogenetic Species Variability (PSV)
## Análise com dados de presença e ausência das espécies nas comunidades.
resultados_PSV <- psv(composicao_especies_P,filogenia_aves)

## Mostra os valores de PSV para cada comunidade.
resultados_PSV 

# Regularidade da diversidade alfa filogenética
# Variance of Pairwise Distance (VPD)
## Transformando data frame em matriz.
dados_matriz <- as.matrix(composicao_especies_P)

## Transformar os dados para o formato requerido pelo pacote pez.
dados <- comparative.comm(filogenia_aves, dados_matriz)

## Análise.
resultados_VPD <- .vpd(dados, cophenetic(filogenia_aves))

## Mostra os valores de VPD para cada comunidade.
resultados_VPD 

# Correlação entre as métricas de diversidade alfa filogenética
## Data frame
# Vamos criar um data.frame com os resultados das métricas da dimensão riqueza.
metricas_riqueza <- data.frame(riqueza = resultados_PD$SR,
                               PD = resultados_PD$PD,
                               PSR = resultados_PSR$PSR,
                               PE = resultados_PE)

## Gráfico
# Gráfico mostrando na parte:
# i) inferior a distribuição dos pontos considerando as métricas pareadas
# ii) superior o valor da correlação de pearson
# iii) diagonal a curva de densidade
ggpairs(metricas_riqueza, upper = list(continuous = wrap("cor", size = 4))) +
  tema_livro()

## Data frame
# Vamos criar um data.frame com os resultados das métricas da dimensão divergência.
metricas_divergencia <- data.frame(riqueza = resultados_PD$SR,
                                   MPD = resultados_MPD_PA,
                                   MPD_AB = resultados_MPD_AB,
                                   MNTD = resultados_MNTD_PA,
                                   MNTD_AB = resultados_MNTD_AB,
                                   PSV = resultados_PSV$PSVs)

## Gráfico
ggpairs(metricas_divergencia, upper = list(continuous = wrap("cor", size = 4))) +
  tema_livro()

# Associação entre a diversidade alfa filogenética e o ambiente
## Dados
# Vamos inserir os dados de precipitação na planilha metrica_divergencia.
metricas_divergencia$precipitacao <- precipitacao_filogenetica$prec

## Gráficos
MPD_PA_plot <- ggplot(metricas_divergencia, aes(precipitacao, MPD)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  labs(x = "Precipitação (mm)", 
       y = "Mean Pairwise Distance\n (MPD - Ausência e Presença)") +
  tema_livro()

MPD_AB_plot <- ggplot(metricas_divergencia, aes(precipitacao, MPD_AB)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  labs(x = "Precipitação (mm)", 
       y = "Mean Pairwise Distance\n (MPD - Abundância)", size = 8) +
  tema_livro() 

MNTD_AP_plot <- ggplot(metricas_divergencia, aes(precipitacao, MNTD)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  labs(x = "Precipitação (mm)", 
       y = "Mean Nearest Taxon Distance\n (MNTD - Ausência e Presença)", 
       size = 8) +
  tema_livro() 

MNTD_AB_plot <- ggplot(metricas_divergencia, aes(precipitacao, MNTD_AB)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Precipitação (mm)", 
       y = "Mean Nearest Taxon Distance\n (MNTD - Abundância)", 
       size = 8) +
  tema_livro()

ggarrange(MPD_PA_plot, MPD_AB_plot, MNTD_AP_plot, MNTD_AB_plot,
          ncol = 2, nrow = 2)

##  gráficos das métricas da dimensão riqueza da diversidade alfa filogenética
## Dados
# Vamos inserir os dados de precipitação na planilha metrica_riqueza.
metricas_riqueza$precipitacao <- precipitacao$prec

## Gráficos
Riqueza_plot <- ggplot(metricas_riqueza, aes(precipitacao, riqueza)) +
  geom_point(size = 4, shape = 19, col = "darkorange") +
  geom_smooth(method = lm, se = FALSE, color = "black") + 
  labs(x = "Precipitação (mm)", y = "Riqueza de espécies") +
  tema_livro() 

PD_plot <- ggplot(metricas_riqueza, aes(precipitacao, PD)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Precipitação (mm)", 
       y = "Diversidade Filogenética\n (Faith)", size = 8) +
  tema_livro()

PSR_plot <- ggplot(metricas_riqueza, aes(precipitacao, PSR)) +
  geom_point(size = 4, shape = 19, col = "darkorange")  +
  geom_smooth(method = lm, se = FALSE, color = "black") + 
  labs(x = "Precipitação (mm)", 
       y = "Phylogenetic Species Richness\n (PSR)", 
       size = 8) +
  tema_livro()

PE_plot <- ggplot(metricas_riqueza, aes(precipitacao, PE)) +
  geom_point(size = 4, shape = 19, col = "darkorange")  +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Precipitação (mm)", 
       y = "Phylogenetic Endemism\n (PE)", 
       size = 8) + 
  tema_livro() 

ggarrange(Riqueza_plot, PD_plot, PSR_plot, PE_plot, ncol = 2, nrow = 2)

# Métricas de diversidade beta filogenética
# Divergência da diversidade beta filogenética
# Community Mean Pairwise Distance (COMDIST)
# Análise com dados de presença e ausência das espécies nas comunidades.
resultados_Comdist_PA <- comdist(composicao_especies_P, 
                                 cophenetic(filogenia_aves), 
                                 abundance.weighted = FALSE)
resultados_Comdist_PA

# Análise com dados de abundância das espécies nas comunidades.
resultados_Comdist_AB <- comdist(composicao_especies_P, 
                                 cophenetic(filogenia_aves), 
                                 abundance.weighted = TRUE)
resultados_Comdist_AB

# Community Mean Nearest Taxon Distance (COMDISTNT)
## Análise com dados de presença e ausência das espécies nas comunidades.
resultados_Comdistnt_PA <- comdistnt(composicao_especies_P, 
                                     cophenetic(filogenia_aves), 
                                     abundance.weighted = FALSE)
resultados_Comdistnt_PA

## Análise com dados de abundância das espécies nas comunidades.
resultados_Comdistnt_AB <- comdistnt(composicao_especies_P, 
                                     cophenetic(filogenia_aves), 
                                     abundance.weighted = TRUE)
resultados_Comdistnt_AB

# Correlação entre as métricas de diversidade beta filogenética

## Dados
## Vamos criar um data frame com os resultados das métricas da dimensão divergência.
metricas_divergencia_beta <- data.frame(
  COMDIST_PA = as.numeric(resultados_Comdist_PA),
  COMDIST_AB = as.numeric(resultados_Comdist_AB),
  COMDISTNT_PA = as.numeric(resultados_Comdistnt_PA),
  COMDISTNT_AB = as.numeric(resultados_Comdistnt_AB))

## Gráfico
ggpairs(metricas_divergencia_beta,
        upper = list(continuous = wrap("cor", size = 4))) +
  tema_livro()

# Associação entre a divergência da diversidade beta filogenética e o ambiente
## Dados
## Precisamos calcular a dissimilaridade par a par da precipitação entre as comunidades.
dis_prec <- vegdist(precipitacao, "euclidian")

## Vamos inserir estes dados na planilha metrica_divergencia_beta.
metricas_divergencia_beta$dis_prec <- as.numeric(dis_prec)

## Gráficos.
COMDIST_PA_plot <- ggplot(metricas_divergencia_beta, 
                          aes(dis_prec, COMDIST_PA)) +
  geom_point(size = 4, shape = 19, col = "darkorange") +
  labs(x = "Diferença na precipitação (mm)", 
       y = "COMDIST\n (Presença e Ausência)") + 
  tema_livro() 

COMDIST_AB_plot <- ggplot(metricas_divergencia_beta, 
                          aes(dis_prec, COMDIST_AB)) +
  geom_point(size = 4, shape = 19, col = "darkorange") +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Diferença na precipitação (mm)", 
       y = "COMDIST\n (Abundância)", size = 8) +
  tema_livro()

COMDISTNT_PA_plot <- ggplot(metricas_divergencia_beta, 
                            aes(dis_prec, COMDISTNT_PA)) +
  geom_point(size = 4, shape = 19, col = "darkorange") +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Diferença na precipitação (mm)", 
       y = "COMDISTNT\n (Ausência e Presença)", 
       size = 8) + 
  tema_livro()

COMDISTNT_AB_plot <- ggplot(metricas_divergencia_beta, 
                            aes(dis_prec, COMDISTNT_AB)) +
  geom_point(size = 4, shape = 19, col = "darkorange") +
  labs(x = "Diferença na precipitação (mm)", 
       y = " COMDISTNT\n (Abundância)", 
       size = 8) +
  tema_livro() 

ggarrange(COMDIST_PA_plot, COMDIST_AB_plot, COMDISTNT_PA_plot, 
          COMDISTNT_AB_plot, ncol = 2, nrow = 2)

# Riqueza da diversidade beta filogenética
# Phylogenetic index of beta diversity (Phylosor)
## Análise com dados de presença e ausência das espécies nas comunidades.
resultados_Phylosor <- picante::phylosor(composicao_especies_P, filogenia_aves)

## Mostra uma matriz triangular com a similaridade entre a fração dos ramos 
## compartilahdos entre duas comunidades
resultados_Phylosor 

# Unique Fraction metric (UniFrac)
## Análise com dados de presença e ausência das espécies nas comunidades.
resultados_UniFrac <- picante::unifrac(composicao_especies_P, filogenia_aves)
resultados_UniFrac

# Correlação entre Phylosor e Unifrac
## Dados
## Vamos criar um data.frame com os resultados das métricas separados 
## para as dimensões de riqueza e divergência.
metricas_riqueza_beta <- data.frame(Phylosor = as.numeric(resultados_Phylosor),
                                    UniFrac = as.numeric(resultados_UniFrac))

## Gráfico
ggpairs(metricas_riqueza_beta, upper=list(continuous = wrap("cor", size = 4))) +
  tema_livro()

# Associação entre a riqueza da diversidade beta filogenética e o ambiente
## Dados
# Vamos inserir os dados de precipitação na planilha metrica_riqueza_beta.
metricas_riqueza_beta$dis_prec <- as.numeric(dis_prec)

## Gráficos
## Phylosor.
plot_phylosor <- ggplot(metricas_riqueza_beta, aes(dis_prec, Phylosor)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  scale_y_continuous(limits = c(0, 1.0)) +
  labs(x = "Diferença na precipitação (mm)", 
       y = "Phylosor", size = 8) +
  tema_livro()

## Unifrag.
plot_unifrac <- ggplot(metricas_riqueza_beta, aes(dis_prec, UniFrac)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  scale_y_continuous(limits = c(0, 1.0)) +
  labs(x = "Diferença na precipitação (mm)", 
       y = "UniFrac", size = 8) +
  tema_livro()

ggarrange(plot_phylosor, plot_unifrac, ncol = 2)

#Partição da diversidade beta filogenética
## Temos que transformar os dados para presença e ausência das espécies nas comunidades.
dados_PA <- decostand(composicao_especies_P, "pa")

## Partição dos componentes do Phylosor.
resultados_Phylosor_particao <- phylo.beta.pair(dados_PA,
                                                filogenia_aves, 
                                                index.family = "sorensen")
## Partição dos componentes do UniFrac.
resultados_UniFrac_particao <- phylo.beta.pair(dados_PA,
                                               filogenia_aves, 
                                               
                                               index.family = "jaccard")
## Gráfico com os resultados dos componentes substituição e aninhamento da diversidade beta filogenética - Phylosor 
## Dados
## Vamos preparar os dados para o gráfico.
particao_phylosor <- data.frame(
  substituicao = as.numeric(resultados_Phylosor_particao$phylo.beta.sim),
  aninhamento = as.numeric(resultados_Phylosor_particao$phylo.beta.sne),
  sorensen = as.numeric(resultados_Phylosor_particao$phylo.beta.sor),
  dis_prec = as.numeric(dis_prec))

## Gráficos
sorensen_plot <- ggplot(particao_phylosor, 
                        aes(dis_prec, sorensen)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "", y = "Sorensen") +
  tema_livro()

subst_plot <- ggplot(particao_phylosor, 
                     aes(dis_prec, substituicao)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Diferença na precipitação\n (mm)", 
       y = "Componente Substituição", size = 8) +
  tema_livro()

aninha_plot <- ggplot(particao_phylosor, 
                      aes(dis_prec, aninhamento)) +
  geom_point(size = 4, shape = 19, col = "darkorange") + 
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "", y = "Componente aninhamento", size = 8) +
  tema_livro()

ggarrange(sorensen_plot, subst_plot, aninha_plot, 
          ncol = 3, nrow = 1)

# Modelos Nulos
## Nearest Relative Index (NRI) ou Standardized Effect Size of MPD
### NRI ou SES_MPD
resultados_SES_MPD <- ses.mpd(composicao_especies_P, cophenetic(filogenia_aves),
                              null.model = "taxa.labels", 
                              abundance.weighted = FALSE,
                              runs = 999) 

### Mostra a riqueza de espéices, MPD observado, média e desvio padrão dos 
### valores de MPD das aleatorizações, SES e o valor de p.
### head(resultados_SES_MPD)

## Nearest Taxon Index (NTI) ou Standardized Effect Size of MNTD
## NTI ou SES_MNTD
resultados_SES_MNTD <- ses.mntd(composicao_especies_P, cophenetic(filogenia_aves),
                                null.model = "taxa.labels", 
                                abundance.weighted = FALSE,
                                runs = 999) 

### Mostra a riqueza de espéices,MNTD observado, média e desvio padrão dos 
### valores de MNTD das aleatorizações, SES e o valor de p.
head(resultados_SES_MNTD)

# Standardized Effect Size of PD
### SES_PD
resultados_SES_PD <- ses.pd(composicao_especies_P, filogenia_aves,
                            null.model = "independentswap", 
                            runs = 999) 

### Mostra a riqueza de espéices,MNTD observado, média e desvio padrão dos 
### valores de PD das aleatorizações, SES e o valor de p.
head(resultados_SES_PD) 

## Standardized effect size do Phylosor
### Modelo nulo que rearranja o nome das espécies na filogenia. 
modelos_nulo <- phylosor.rnd(composicao_especies_P, filogenia_aves, 
                             null.model = "taxa.labels", runs = 9)

### Função para calcular o SES eo valor de P.
ses.physo <- function(obs, nulo_phylosor){
  nulo_phylosor <- t(as.data.frame(lapply
                                   (nulo_phylosor, as.vector)))
  physo.obs <- as.numeric(obs)
  physo.mean <- apply(nulo_phylosor, MARGIN = 2, 
                      FUN = mean, na.rm = TRUE)
  physo.sd <- apply(nulo_phylosor, MARGIN = 2, 
                    FUN = sd, na.rm = TRUE)
  physo.ses <- (physo.obs - physo.mean)/physo.sd
  physo.obs.rank <- apply(X = rbind(physo.obs, 
                                    nulo_phylosor), MARGIN = 2, 
                          FUN = rank)[1, ]
  physo.obs.rank <- ifelse(is.na(physo.mean), NA, 
                           physo.obs.rank)
  data.frame(physo.obs, physo.mean, physo.sd, 
             physo.obs.rank, physo.ses, 
             physo.obs.p = physo.obs.rank/
               (dim(nulo_phylosor)[1] + 1))
}

### Resultados
resultados <- ses.physo (resultados_Phylosor, modelos_nulo)
head(resultados)

