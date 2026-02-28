## Pacotes
library(ecodados)
library(here)
library(tidyverse)
library(sf)
library(raster)
library(rgdal)
library(spData)
library(rnaturalearth)
library(geobr)
library(ggplot2)
library(ggspatial)
library(tmap)
library(tmaptools)
library(grid)
library(mapview)
library(leaflet)
library(viridis)
library(knitr)
library(sidrar)
library(landscapetools)
library(colorspace)

## Dados
world <- world
volcano <- volcano
geo_anfibios_locais <- ecodados::geo_anfibios_locais
geo_anfibios_especies <- ecodados::geo_anfibios_especies
geo_vetor_nascentes <- ecodados::geo_vetor_nascentes
geo_vetor_hidrografia <- ecodados::geo_vetor_hidrografia
geo_vetor_cobertura <- ecodados::geo_vetor_cobertura
geo_vetor_rio_claro <- ecodados::geo_vetor_rio_claro
geo_vetor_brasil <- ecodados::geo_vetor_brasil
geo_vetor_brasil_anos <- ecodados::geo_vetor_brasil_anos
geo_vetor_am_sul <- ecodados::geo_vetor_am_sul
geo_vetor_biomas <- ecodados::geo_vetor_biomas
geo_vetor_mata_atlantica <- ecodados::geo_vetor_mata_atlantica
geo_raster_srtm <- ecodados::geo_raster_srtm
geo_raster_bioclim <- ecodados::geo_raster_bioclim
geo_raster_globcover_mata_atlantica <- ecodados::geo_raster_globcover_mata_atlantica


# Vetor

## sf: principal pacote no R para dados vetoriais
## Dados vetoriais de polígonos do mundo
data(world)
world

## Plot dos polígonos do mundo
plot(world[1], col = viridis::viridis(100), main = "Mapa do mundo")

# Raster
##  raster: principal pacote no R para dados raster
## Dados de altitude de um vulcão
volcano[1:5, 1:5]

## Vamos transformar essa matriz de dados em um raster com a função raster::raster().
## Rasterlayer
raster_layer <- raster::raster(volcano)
raster_layer

## Plot raster layers
plot(raster_layer, col = viridis::viridis(n = 100))

## Raster layers
raster_layer1 <- raster_layer
raster_layer2 <- raster_layer * raster_layer
raster_layer3 <- sqrt(raster_layer)
raster_layer4 <- log10(raster_layer)

## Raster brick
raster_brick <- raster::brick(
    raster_layer1, raster_layer2,
    raster_layer3, raster_layer4
)
raster_brick

## Plot raster brick
plot(raster_brick, col = viridis::viridis(n = 25), main = "")

## Raster layers
raster_layer1 <- raster_layer
raster_layer2 <- raster_layer * raster_layer
raster_layer3 <- sqrt(raster_layer)
raster_layer4 <- log10(raster_layer)

## Raster stack
raster_stack <- raster::stack(
    raster_layer1, raster_layer2,
    raster_layer3, raster_layer4
)
raster_stack
## Plot raster stack
plot(raster_stack, col = viridis::viridis(n = 25), main = "")

# Sistema de Referência de Coordenadas e Unidades
## Sistema de Referência de Coordenadas (CRS) no R
## Listagem dos Sistemas de Referências de Coordenadas no R
crs_data <- rgdal::make_EPSG()
head(crs_data)

# Importar dados
## Formatos vetoriais importados e exportados pelo pacote sf
head(sf::st_drivers())
## Criar diretório
dir.create(here::here("dados"))
dir.create(here::here("dados", "vetor"))
## Aumentar o tempo de download
options(timeout = 1e3)

## Download
for (i in c(".dbf", ".prj", ".shp", ".shx")) {
    # Pontos de nascentes
    download.file(
        url = paste0("http://geo.fbds.org.br/SP/RIO_CLARO/HIDROGRAFIA/SP_3543907_NASCENTES", i),
        destfile = here::here("dados", "vetor", paste0("SP_3543907_NASCENTES", i)), mode = "wb"
    )

    # Linhas de hidrografia
    download.file(
        url = paste0("http://geo.fbds.org.br/SP/RIO_CLARO/HIDROGRAFIA/SP_3543907_RIOS_SIMPLES", i),
        destfile = here::here("dados", "vetor", paste0("SP_3543907_RIOS_SIMPLES", i)), mode = "wb"
    )

    # Polígonos de cobertura da terra
    download.file(
        url = paste0("http://geo.fbds.org.br/SP/RIO_CLARO/USO/SP_3543907_USO", i),
        destfile = here::here("dados", "vetor", paste0("SP_3543907_USO", i)), mode = "wb"
    )
}

## Importar nascentes
geo_vetor_nascentes <- sf::st_read(
    here::here("dados", "vetor", "SP_3543907_NASCENTES.shp"),
    quiet = TRUE
)
## Plot
plot(geo_vetor_nascentes[1],
    pch = 20, col = "blue", main = NA,
    axes = TRUE, graticule = TRUE
)

## Importar hidrografia
geo_vetor_hidrografia <- sf::st_read(
    here::here("dados", "vetor", "SP_3543907_RIOS_SIMPLES.shp"),
    quiet = TRUE
)
## Plot
plot(geo_vetor_hidrografia[1], col = "steelblue", main = NA, axes = TRUE, graticule = TRUE)

## Importar cobertura da terra
geo_vetor_cobertura <- sf::st_read(
    here::here("dados", "vetor", "SP_3543907_USO.shp"),
    quiet = TRUE
)
## Plot
plot(geo_vetor_cobertura[5],
    pal = colorspace::qualitative_hcl,
    main = NA, axes = TRUE, graticule = TRUE,
    key.pos = NULL, # Remove a legenda lateral padrão que usa o 'layout'
    reset = FALSE # Evita erros de reset do painel depois
)
legend(
    "left",
    inset = 0.05, pch = 15, cex = .7, pt.cex = 2.5,
    legend = levels(factor(geo_vetor_cobertura$CLASSE_USO)),
    col = colorspace::qualitative_hcl(length(levels(factor(geo_vetor_cobertura$CLASSE_USO))))
)

## Importar utilizando pacotes
## Listar todos os dados do geobr
geobr::list_geobr()

## Polígono do limite do município de Rio Claro
geo_vetor_rio_claro <- geobr::read_municipality(code_muni = 3543907, 
                                                year = 2020, showProgress = FALSE)
## Plot
plot(geo_vetor_rio_claro[1], col = "gray", main = NA, axes = TRUE, graticule = TRUE)
## Polígono do limite do Brasil
geo_vetor_brasil <- rnaturalearth::ne_countries(scale = "large", 
                                                country = "Brazil", returnclass = "sf")
## Plot
plot(geo_vetor_brasil[1], col = "gray", main = NA, axes = TRUE, graticule = TRUE)

## Criar diretório
dir.create(here::here("dados", "tabelas"))
## Download
download.file(url = "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2392&file=ecy2392-sup-0001-DataS1.zip",
              destfile = here::here("dados", "tabelas", "atlantic_amphibians.zip"), mode = "wb")

## Unzip
unzip(zipfile = here::here("dados", "tabelas", "atlantic_amphibians.zip"),
      exdir = here::here("dados", "tabelas"))

## Importar os dados pelo pacote ecodados
geo_anfibios_locais <- ecodados::geo_anfibios_locais


## Converter dados tabulares para sf
geo_anfibios_locais_vetor <- geo_anfibios_locais |>
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
geo_anfibios_locais_vetor

## Plot
plot(geo_anfibios_locais_vetor[1], pch = 20, col = "black", 
     main = NA, axes = TRUE, graticule = TRUE)

# Converter dados espaciais sp para sf
## Polígonos países sp
co110_sp <- rnaturalearth::countries110
class(co110_sp)
## Polígonos países sf
co110_sf <- sf::st_as_sf(co110_sp)
class(co110_sf)
## Polígonos países sp
co110_sp <- sf::as_Spatial(co110_sf)
class(co110_sp)

## Criar diretório
dir.create(here::here("dados", "raster"))

## Aumentar o tempo de download
options(timeout = 1e3)

## Download
download.file(url = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_27_17.zip",
              destfile = here::here("dados", "raster", "srtm_27_17.zip"), mode = "wb")

## Unzip
unzip(zipfile = here::here("dados", "raster", "srtm_27_17.zip"),
      exdir = here::here("dados", "raster"))

## Importar raster de altitude
geo_raster_srtm <- raster::raster(here::here("dados", "raster", "srtm_27_17.tif"))
geo_raster_srtm

