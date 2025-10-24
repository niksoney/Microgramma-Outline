# =============================================================================
# PGLS
# =============================================================================
# Article: ........................... 
# Authors: Niksoney A. Mendonça, Juliana Aljahara, D. Victor S. Silva and Thaís E. Almeida
# Script author: Niksoney Azevedo Mendonça
# E-mail: niksoneyazevedo2017@gmail.com
# =============================================================================
# Instalar e Carregar Pacotes
# =============================================================================
# Instalar todos os pacotes de uma vez
install.packages(c("ape","caper","nlme","phytools","ggtree","patchwork","dplyr","ggplot2","ggnewscale"))

# Carregar todos os pacotes
library(ape)         # Análises e manipulação de árvores filogenéticas
library(caper)       # Modelos comparativos filogenéticos (PGLS, independência)
library(nlme)        # Modelos lineares e não lineares mistos
library(phytools)    # Ferramentas adicionais para análises filogenéticas
library(ggtree)      # Visualização de árvores com ggplot2
library(patchwork)   # Combinação e layout de múltiplos gráficos ggplot2
library(dplyr)       # Manipulação eficiente de dados
library(ggplot2)     # Gráficos e visualizações
library(ggnewscale)  # Adicionar múltiplas escalas de cores no ggplot2

#Definir meu diretório
setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2.Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine/phylogenetic generalized least squares regressions")

sterile_spectrum <- read.csv("espectro_esteril.csv")
espectro_fertil <- read.csv("espectro_fertil.csv")
espectros <- read.csv("espectros.csv")
contorno_esteril <- read.csv("contorno_esteril.csv")
contorno_fertil <- read.csv("contorno_fertil.csv")
contornos <- read.csv("contornos.csv")

# Substituir nomes na tabela contorno_esteril
contorno_esteril$especie <- recode(contorno_esteril$especie,
                                   "DIC" = "Microgramma dictyophylla", 
                                   "LAT" = "Microgramma latevagans", 
                                   "NAN" = "Microgramma nana", 
                                   "PIL" = "Microgramma piloselloides", 
                                   "TOB" = "Microgramma tobagensis", 
                                   "TEC" = "Microgramma tecta", 
                                   "REP" = "Microgramma reptans", 
                                   "PER" = "Microgramma percussa"
)

# Substituir nomes na tabela contorno_fertil
contorno_fertil$especie <- recode(contorno_fertil$especie,
                                  "DIC" = "Microgramma dictyophylla", 
                                  "LAT" = "Microgramma latevagans", 
                                  "NAN" = "Microgramma nana", 
                                  "PIL" = "Microgramma piloselloides", 
                                  "TOB" = "Microgramma tobagensis", 
                                  "TEC" = "Microgramma tecta", 
                                  "REP" = "Microgramma reptans", 
                                  "PER" = "Microgramma percussa"
)

# Substituir nomes na tabela contornos
contornos$especie <- recode(contornos$especie,
                            "DIC" = "Microgramma dictyophylla", 
                            "LAT" = "Microgramma latevagans", 
                            "NAN" = "Microgramma nana", 
                            "PIL" = "Microgramma piloselloides", 
                            "TOB" = "Microgramma tobagensis", 
                            "TEC" = "Microgramma tecta", 
                            "REP" = "Microgramma reptans", 
                            "PER" = "Microgramma percussa"
)


# Função ajustada: roda PCA, mostra resultado e salva os PCs que explicam ≥ 95%
ver_pca_95_salvar <- function(tabela, nome, arquivo_saida) {
  dados_num <- tabela[,-1]  # Remove a coluna com os nomes das espécies
  pca <- prcomp(dados_num, scale. = TRUE)
  
  # Variância acumulada
  variancia_acum <- summary(pca)$importance[3,]
  
  # Quantos PCs explicam ≥ 95%
  num_pcs_95 <- which(variancia_acum >= 0.95)[1]
  
  # Mensagem no console
  cat(paste0("\n📊 ", nome, ": ", num_pcs_95, " PCs explicam ≥ 95% da variância.\n"))
  cat("Variância acumulada dos primeiros 10 PCs:\n")
  print(round(variancia_acum[1:10], 3))
  
  # Criar dataframe com os escores desejados
  escores_filtrados <- data.frame(
    especie = tabela[,1],
    pca$x[, 1:num_pcs_95]
  )
  
  # Salvar em CSV
  write.csv(escores_filtrados, arquivo_saida, row.names = FALSE)
  
  return(pca)
}

pca_espectro_esteril <- ver_pca_95_salvar(espectro_esteril, "Espectro Estéril", "pca_espectro_esteril_90.csv")
pca_espectro_fertil  <- ver_pca_95_salvar(espectro_fertil,  "Espectro Fértil",  "pca_espectro_fertil_90.csv")
pca_espectros        <- ver_pca_95_salvar(espectros,        "Espectros (Todos)", "pca_espectros_90.csv")
pca_contorno_esteril <- ver_pca_95_salvar(contorno_esteril, "Contorno Estéril", "pca_contorno_esteril_90.csv")
pca_contorno_fertil  <- ver_pca_95_salvar(contorno_fertil,  "Contorno Fértil",  "pca_contorno_fertil_90.csv")
pca_contornos        <- ver_pca_95_salvar(contornos,        "Contornos (Todos)", "pca_contornos_90.csv")

##########################################################################################################
##########################################################################################################
##########################Phylogenetic generalized least squares regressions##############################
##########################################################################################################
##########################################################################################################

rm(list = ls())  # Remove todos os objetos do ambiente global 

# Definir diretório
setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2.Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine/PGLs")

# Carregar árvore
arvore <- read.tree("scaly_podada - Copia.newick")
plot(arvore, type = "phylogram", show.tip.label = TRUE)

# Lista com nomes dos arquivos e número de PCs a usar
datasets_info <- list(
  sterile_spectrum = list(file = "pca_espectro_esteril_90.csv", n_pcs = 2),
  sterile_shape = list(file = "pca_contorno_esteril_90.csv", n_pcs = 2),
  fertile_spectrum = list(file = "pca_espectro_fertil_90.csv",  n_pcs = 2),
  fertile_shape  = list(file = "pca_contorno_fertil_90.csv",  n_pcs = 2),
  spectra = list(file = "pca_espectros_90.csv",         n_pcs = 2),
  shape = list(file = "pca_contornos_90.csv",        n_pcs = 2)
)

run_pgls_with_data <- function(data_file, n_pcs, tree) {
  # Carregar dados
  dados <- read.csv(data_file)
  
  # Corrigir nomes: substituir espaços por underscores para compatibilidade com a árvore
  dados$especie <- gsub(" ", "_", dados$especie)
  
  # Verificar e remover possíveis duplicatas
  if(any(duplicated(dados$especie))) {
    warning("Espécies duplicadas encontradas. Mantendo apenas a primeira ocorrência.")
    dados <- dados[!duplicated(dados$especie), ]
  }
  
  # Filtrar espécies presentes na árvore e ordenar conforme a árvore
  especies_comuns <- intersect(dados$especie, tree$tip.label)
  if(length(especies_comuns) == 0) {
    stop("Nenhuma espécie em comum entre os dados e a árvore.")
  }
  
  dados <- dados[dados$especie %in% especies_comuns, ]
  dados <- dados[match(tree$tip.label, dados$especie), ]
  
  # Criar objeto comparative.data
  comp <- comparative.data(phy = tree, 
                           data = dados, 
                           names.col = "especie", 
                           vcv = TRUE,
                           warn.dropped = TRUE)
  
  # Criar fórmula
  vars <- paste0("PC", 2:n_pcs)
  formula <- as.formula(paste("PC1 ~", paste(vars, collapse = " + ")))
  
  # Ajustar modelo PGLS
  modelo <- pgls(formula, data = comp, lambda = "ML")
  
  # Retornar resultados
  list(
    modelo = summary(modelo),
    dados = dados,
    comp = comp
  )
}

# Executar análises para todos os datasets
resultados_completos <- lapply(datasets_info, function(info) {
  tryCatch({
    run_pgls_with_data(info$file, info$n_pcs, arvore)
  }, error = function(e) {
    message(paste("Erro ao processar", info$file, ":", e$message))
    return(NULL)
  })
})

# Criar o objeto resultados_completos
resultados_completos <- lapply(datasets_info, function(info) {
  run_pgls_with_data(info$file, info$n_pcs, arvore)
})

# Imprimir os sumários dos modelos PGLS para cada dataset
for (nome in names(resultados_completos)) {
  cat("\n===========================\n")
  cat("Resultado para:", nome, "\n")
  cat("===========================\n")
  print(resultados_completos[[nome]]$modelo)
}

############################################################################
#                                                                          #
#      GRÁFICO DE RESULTADOS PGLS                                          #
#                                                                          #
############################################################################

# --- Passo 1: Carregar Pacotes e Definir Cores ---
# Certifique-se de que esses pacotes estão instalados e carregados
library(ggplot2)
library(dplyr)
library(patchwork) # Para combinar os gráficos

# Paleta de cores que já definimos
nature_palette <- c(
  "forest" = "#1A5225",
  "sea" = "#2A7A7E",
  "moss" = "#6A8A28",
  "gold_significant" = "#D4AF37", # Dourado para regressão significativa
  "gray_nonsignificant" = "#888888"  # Cinza para não significativa
)

# --- VERSÃO FINAL: CONSTELAÇÃO COM LINHAS COLORIDAS POR ESPÉCIE ---

plot_pgls_regression <- function(resultado, titulo) {
  
  # Extrair dados e modelo
  dados <- resultado$dados
  modelo_summary <- resultado$modelo
  
  # Extrair coeficientes e estatísticas
  intercepto <- modelo_summary$coefficients["(Intercept)", "Estimate"]
  inclinacao <- modelo_summary$coefficients["PC2", "Estimate"]
  p_valor <- modelo_summary$coefficients["PC2", "Pr(>|t|)"]
  r_squared <- round(modelo_summary$r.squared, 3)
  lambda <- round(modelo_summary$param["lambda"], 3)
  
  # Formatar p-valor
  p_texto <- ifelse(p_valor < 0.001, "< 0.001", as.character(round(p_valor, 3)))
  
  # Definir cor da linha de regressão
  cor_linha <- ifelse(p_valor < 0.05, nature_palette["gold_significant"], nature_palette["gray_nonsignificant"])
  
  # Criar a expressão para o subtítulo
  subtitle_expression <- bquote(
    "R"^2 == .(r_squared) ~ " | " ~ italic(P) == .(p_texto) ~ " | " ~ lambda == .(lambda)
  )
  
  # Calcular o centroide (ponto médio) do morfoespaço
  centroide_x <- mean(dados$PC2, na.rm = TRUE)
  centroide_y <- mean(dados$PC1, na.rm = TRUE)
  
  # --- GRÁFICO COM A NOVA ESTÉTICA ---
  p <- ggplot(dados, aes(x = PC2, y = PC1)) +
    
    # EFEITO CONSTELAÇÃO: Linhas de conexão com a cor da espécie
    geom_segment(
      aes(xend = centroide_x, yend = centroide_y, color = gsub("_", " ", especie)), # MUDANÇA AQUI
      alpha = 0.6,      # Opacidade para as linhas não dominarem
      linewidth = 0.8,
      show.legend = FALSE # Importante para não criar uma legenda duplicada
    ) +
    
    # Pontos principais (por cima das linhas)
    geom_point(
      aes(color = gsub("_", " ", especie)),
      size = 12
    ) +
    
    # Linha de regressão PGLS
    geom_abline(
      intercept = intercepto,
      slope = inclinacao,
      color = cor_linha,
      linewidth = 1.5
    ) +
    
    # Títulos e o subtítulo estatístico
    labs(
      title = gsub("_", " ", titulo),
      subtitle = subtitle_expression,
      x = "PC2",
      y = "PC1", face = "bold",
      color = NULL # Remove o título da legenda
    ) +
    
    # Tema claro e profissional
    theme_classic(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5, color = "gray30"),
      axis.title = element_text(size = 22),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray95", linetype = "dotted"),
      plot.background = element_rect(fill = "#FFFFFF", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA)
    )
  
  return(p)
}

# --- Passo 3: Gerar e Combinar os Gráficos ---

# Paleta de cores para as espécies (certifique-se que os nomes estão corretos)
cores_especies <- c(
  "Microgramma tecta" = "#CC79A7", "Microgramma nana" = "#009E73",
  "Microgramma piloselloides" = "#0072B2", "Microgramma percussa" = "#F0E442",
  "Microgramma tobagensis" = "#555555", "Microgramma latevagans" = "#56B4E9",
  "Microgramma dictyophylla" = "#E69F00", "Microgramma reptans" = "#D55E00"
)

# Criar uma lista para armazenar os gráficos
lista_de_plots <- list()

# Loop através dos resultados PGLS para criar cada gráfico
for (nome in names(resultados_completos)) {
  if (!is.null(resultados_completos[[nome]])) {
    p <- plot_pgls_regression(resultados_completos[[nome]], nome)
    lista_de_plots[[nome]] <- p
  }
}

# --- COMBINAR GRÁFICOS E AJUSTAR LEGENDA FINAL PARA TEMA CLARO ---

# (O loop que cria a 'lista_de_plots' continua o mesmo)

# Combinar todos os gráficos
grafico_final_pgls <- wrap_plots(lista_de_plots, ncol = 2, nrow = 3) +
  plot_layout(guides = "collect") &
  scale_color_manual(values = cores_especies) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "italic", size = 18)
  )

# Exibir o gráfico final
print(grafico_final_pgls)

# Salvar o gráfico final
ggsave("Imagens_prontas/PGLS_Final_Panel.png", grafico_final_pgls, width = 16, height = 18, dpi = 300)
