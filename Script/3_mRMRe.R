# =============================================================================
# mRMRe
# =============================================================================
# Article: ........................... 
# Authors: Niksoney A. Mendonça, Juliana Aljahara, D. Victor S. Silva e Thaís E. Almeida
# Script author: Niksoney Azevedo Mendonça
# E-mail: niksoneyazevedo2017@gmail.com
# =============================================================================
# Instalar e Carregar Pacotes
# =============================================================================
# Instalar todos os pacotes de uma vez
install.packages(c("ggplot2","dplyr","mRMRe"))

# Carregar todos os pacotes
library(ggplot2)   # Criação de gráficos elegantes
library(dplyr)     # Manipulação de dados com sintaxe intuitiva
library(mRMRe)     # Seleção e priorização de variáveis por informação mútua

#######################################################################################################
###################################SEPARAR E CONTRUIR AS MATRIZES######################################
#######################################################################################################

#Diretório 
setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2.Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine")

# 1. Carregar as matrizes 
Matriz <- read.csv("Matriz_NIR.csv", sep = ";")
colunas_remover <- c(1, 3, 4, 6, 7, 8, 10, 11)
Matriz <- Matriz[, -colunas_remover]

# 3. Substituir vírgula por ponto em todas as colunas numéricas
Matriz[, 4:ncol(Matriz)] <- lapply(Matriz[, 4:ncol(Matriz)], function(x) {
  if (is.numeric(x) || is.character(x)) {
    as.numeric(gsub(",", ".", as.character(x)))
  } else {
    x
  }
})

# Dividir a matriz em Fertile e Estéril com base na coluna 2
# Assumindo que a coluna 2 contém "fertile" e "sterile" ou similar
Matriz_Fertile <- Matriz[Matriz[, 2] == "Fertile leaf", ]
Matriz_Sterile <- Matriz[Matriz[, 2] == "Sterile leaf", ]

# Verificar os resultados
head(Matriz_Fertile)
head(Matriz_Sterile)

# Contar o número de amostras em cada grupo
nrow(Matriz_Fertile)
nrow(Matriz_Sterile)

# Salvar como CSV (com ponto e vírgula como separador, por exemplo)
write.csv2(Matriz_Fertile, file = "Matriz_Fertile.csv", row.names = FALSE)  # Usa "," para decimais e ";" 
write.csv2(Matriz_Sterile, file = "Matriz_Sterile.csv", row.names = FALSE)  # Usa "," para decimais e ";" 

#######################################################################################################
#FAZER A ANÁLISE mRMRe#################################################################################
#######################################################################################################

rm(list = ls())  # Remove todos os objetos do ambiente global  

setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2.Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine")

Matrizfer <- read.csv("Matriz_Fertile.csv", sep = ";")
Matrizste <- read.csv("Matriz_Sterile.csv", sep = ";")
Matrizcomb <- read.csv("Matriz_NIR.csv", sep = ";")

# Função para substituir vírgulas por pontos
substituir_virgulas <- function(x) {
  if (is.character(x)) {
    as.numeric(gsub(",", ".", x))
  } else {
    x
  }
}

# Função para encontrar o cotovelo geométrico
find_elbow <- function(y) {
  n <- length(y)
  first <- c(1, y[1])
  last <- c(n, y[n])
  line_vec <- last - first
  
  distances <- sapply(1:n, function(i) {
    point <- c(i, y[i])
    vec <- point - first
    d <- abs((line_vec[2] * vec[1] - line_vec[1] * vec[2])) / sqrt(sum(line_vec^2))
    return(d)
  })
  
  which.max(distances)
}

# Listar os datasets e nomes
matrizes <- list(
  Fertile = Matrizfer,
  Sterile = Matrizste,
  Combined = Matrizcomb
)

# Lista para guardar resultados
resultados <- list()
cores <- c("blue", "red", "black")  # cores diferentes para o gráfico
linetypes <- c("solid", "dashed", "dotted")

# Loop por cada matriz
for (i in seq_along(matrizes)) {
  nome <- names(matrizes)[i]
  mat <- matrizes[[i]]
  
  # Excluir colunas 2 e 3, substituir vírgulas
  mat <- mat[, -c(2, 3)]
  mat[, 2:ncol(mat)] <- lapply(mat[, 2:ncol(mat)], substituir_virgulas)
  
  # Calcular média por espécie
  medias <- mat %>%
    group_by(Species) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
  
  # Separar variáveis
  Y <- as.factor(medias$Species)
  X <- medias[, -1]
  
  # Preparar dados para mRMRe
  data_for_mrmr <- data.frame(Y = as.numeric(Y), X)
  dd <- mRMR.data(data = data_for_mrmr)
  
  # Executar mRMR
  set.seed(123)
  feature_selection <- mRMR.classic(data = dd, target_indices = 1, feature_count = 10)
  
  # Calcular R² de todas as bandas
  r2_all <- sapply(colnames(X), function(band) {
    model <- lm(data_for_mrmr[[band]] ~ data_for_mrmr$Y)
    summary(model)$r.squared
  })
  
  ranking <- data.frame(Band = colnames(X), R2 = r2_all)
  ranking <- ranking[order(-ranking$R2), ]
  ranking$Rank <- seq_len(nrow(ranking))
  
  # Encontrar cotovelo
  elbow_index <- find_elbow(ranking$R2)
  bandas_cotovelo <- ranking$Band[1:elbow_index]
  nova_matriz <- medias %>% select(Species, all_of(bandas_cotovelo))
  
  # Armazenar resultados
  resultados[[nome]] <- list(
    ranking = ranking,
    elbow = elbow_index,
    bandas = bandas_cotovelo,
    matriz_final = nova_matriz
  )
}

# --------- Criar gráfico comparativo ---------
# Juntar todos os rankings num único dataframe com coluna 'Matriz'
ranking_completo <- do.call(rbind, lapply(names(resultados), function(nome) {
  df <- resultados[[nome]]$ranking
  df$Matriz <- nome
  df
}))

# Definir cores e linetypes consistentes para cada matriz
matriz_estilos <- data.frame(
  Matriz = names(matrizes),
  Cor = c("blue", "red", "black"),
  Linetype = c("solid", "dashed", "dotted"),
  stringsAsFactors = FALSE
)

# Mesclar com os dados
ranking_completo <- merge(ranking_completo, matriz_estilos, by = "Matriz")

# Criar gráfico base com estilos consistentes
# Criar gráfico base com estilos consistentes
gg <- ggplot(ranking_completo, aes(x = Rank, y = R2, group = Matriz)) +
  geom_line(aes(color = Matriz, linetype = Matriz), size = 1) +
  scale_color_manual(values = setNames(matriz_estilos$Cor, matriz_estilos$Matriz)) +
  scale_linetype_manual(values = setNames(matriz_estilos$Linetype, matriz_estilos$Matriz)) +
  labs(
    title = "R² Comparison by Band Ranking",
    x = "Band Ranking",
    y = "R²",
    color = "Matrices",  # Título da legenda para cor
    linetype = "Matrices"  # Título da legenda para tipo de linha
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )

# Adicionar linhas verticais de cotovelo com estilos consistentes
for (i in seq_along(resultados)) {
  nome <- names(resultados)[i]
  estilo <- matriz_estilos[matriz_estilos$Matriz == nome, ]
  dados <- resultados[[nome]]$ranking
  elbow <- resultados[[nome]]$elbow
  r2_valor <- dados$R2[elbow]
  
  gg <- gg +
    geom_vline(xintercept = elbow, 
               color = estilo$Cor, 
               linetype = estilo$Linetype, 
               size = 1) +
    
    annotate("text",
             x = elbow + 5,
             y = max(dados$R2) * 0.95,
             label = paste0(nome, "\n", elbow, " Spectral bands"),
             color = estilo$Cor,
             size = 5,
             hjust = 0)
}

# Mostrar gráfico
print(gg)

# Salvar gráfico
ggsave(filename = "Imagens_prontas/Comparacao_R2_cotovelo.png",
       plot = gg,
       width = 12, height = 8, units = "in", dpi = 300)

# --------- Salvar as três matrizes finais ---------
write.csv(resultados$Fertile$matriz_final, "Matriz_Fertile_Selecionada.csv", row.names = FALSE)
write.csv(resultados$Sterile$matriz_final, "Matriz_Sterile_Selecionada.csv", row.names = FALSE)
write.csv(resultados$Combined$matriz_final, "Matriz_Combined_Selecionada.csv", row.names = FALSE)

#######################################################################################################
#CONSTRUIR MATRIZES DE DISTÂNCIA#######################################################################
#######################################################################################################

rm(list = ls())  # Remove todos os objetos do ambiente global 

setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2.Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine")

Matrizfer <- read.csv("Matriz_Fertile_Selecionada_ESPECTRO.csv", sep = ",")
Matrizste <- read.csv("Matriz_Sterile_Selecionada_ESPECTRO.csv", sep = ",")
Matrizcomb <- read.csv("Matriz_Combined_Selecionada_ESPECTRO.csv", sep = ",")

# Carregar o pacote dplyr
library(dplyr)

# Definir o mapeamento de nomes
nome_completo <- c(
  "M. dictyophylla" = "Microgramma dictyophylla",
  "M. latevagans" = "Microgramma latevagans",
  "M. nana" = "Microgramma nana",
  "M. piloselloides" = "Microgramma piloselloides",
  "M. tobagensis" = "Microgramma tobagensis",
  "M. tecta" = "Microgramma tecta",
  "M. reptans" = "Microgramma reptans",
  "M. percussa" = "Microgramma percussa"
)

# Aplicar a substituição em cada matriz
Matrizfer <- Matrizfer %>%
  mutate(across(1, ~recode(., !!!nome_completo)))

Matrizste <- Matrizste %>%
  mutate(across(1, ~recode(., !!!nome_completo)))

Matrizcomb <- Matrizcomb %>%
  mutate(across(1, ~recode(., !!!nome_completo)))

# --------- Salvar as três matrizes finais ---------
write.csv(Matrizfer, "espectro_fertil.csv", row.names = FALSE)
write.csv(Matrizste, "espectro_esteril.csv", row.names = FALSE)
write.csv(Matrizcomb, "espectros.csv", row.names = FALSE)
