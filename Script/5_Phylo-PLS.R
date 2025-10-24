# =============================================================================
# Phylo-PLS
# =============================================================================
# Article: ........................... 
# Authors: Niksoney A. Mendonça, Juliana Aljahara, D. Victor S. Silva and Thaís E. Almeida
# Script author: Niksoney Azevedo Mendonça
# E-mail: niksoneyazevedo2017@gmail.com
# =============================================================================
# Instalar e Carregar Pacotes
# =============================================================================
# Instalar todos os pacotes de uma vez
install.packages(c("ape","phytools","dplyr","openxlsx","pls","ggplot2","patchwork","ggphylomorpho"))

# Carregar todos os pacotes
library(ape)           # Funções para árvores filogenéticas
library(phytools)      # Funções adicionais para filogenia e evolução
library(dplyr)         # Manipulação de dados com sintaxe clara
library(openxlsx)      # Leitura/escrita de planilhas Excel
library(pls)           # Análises de regressão PLS
library(ggplot2)       # Criação de gráficos sofisticados
library(patchwork)     # Combinar múltiplos gráficos ggplot2
library(ggphylomorpho) # Funções para representar phylomorphospaces


# 2. Carregar e Preparar os Dados
# =====================================
cat("Carregando todos os conjuntos de dados...\n")
setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2_Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine")

# Função para carregar e limpar dados
load_and_clean_data <- function(filename) {
  df <- read.csv(filename)
  df[, 1] <- gsub("_", " ", df[, 1])
  df[, 1] <- trimws(df[, 1])
  rownames(df) <- df[, 1]
  return(df[, -1, drop = FALSE])
}

# Carregar dados
espectros_geral <- load_and_clean_data("Phylo-PLS/pca_espectros_90.csv")
contornos_geral <- load_and_clean_data("Phylo-PLS/pca_contornos_90.csv")
espectro_esteril <- load_and_clean_data("Phylo-PLS/pca_espectro_esteril_90.csv")
contorno_esteril <- load_and_clean_data("Phylo-PLS/pca_contorno_esteril_90.csv")
espectro_fertil <- load_and_clean_data("Phylo-PLS/pca_espectro_fertil_90.csv")
contorno_fertil <- load_and_clean_data("Phylo-PLS/pca_contorno_fertil_90.csv")

# Carregar árvore filogenética
arvore <- read.tree("Phylo-PLS/scaly_podada - Copia.newick")
arvore$tip.label <- gsub("_", " ", arvore$tip.label)
arvore$tip.label <- trimws(arvore$tip.label)

cat("Dados carregados e preparados.\n")

# =================================================================
# 3. FUNÇÃO DE ANÁLISE (MODIFICADA PARA RETORNAR OS SCORES DO PLS)
# =================================================================
run_phylo_pls_2componentes <- function(dados_forma, dados_espectro, arvore_filogenetica, nome_analise) {
  
  cat(paste("\n--- ANÁLISE:", nome_analise, "---\n"))
  
  # Preparação dos dados (como antes)
  common_species <- intersect(arvore_filogenetica$tip.label, intersect(rownames(dados_forma), rownames(dados_espectro)))
  cat(paste("Espécies em comum:", length(common_species), "\n"))
  if (length(common_species) < 3) { cat("AVISO: Poucas espécies.\n"); return(NULL) }
  
  arvore_podada <- keep.tip(arvore_filogenetica, common_species)
  forma_final <- as.matrix(dados_forma[common_species, c("PC1", "PC2")])
  espectro_final <- as.matrix(dados_espectro[common_species, c("PC1", "PC2")])
  
  forma_padronizada <- scale(forma_final)
  espectro_padronizado <- scale(espectro_final)
  
  C <- vcv.phylo(arvore_podada)[common_species, common_species]
  if (det(C) == 0) { C <- C + diag(nrow(C)) * 1e-6 }
  
  sqrt_invC <- chol(solve(C))
  forma_contrasts <- sqrt_invC %*% forma_padronizada
  espectro_contrasts <- sqrt_invC %*% espectro_padronizado
  
  pls_df <- data.frame(forma_PC1 = forma_contrasts[, 1], forma_PC2 = forma_contrasts[, 2],
                       espectro_PC1 = espectro_contrasts[, 1], espectro_PC2 = espectro_contrasts[, 2])
  
  # Executar PLS
  pls_result <- plsr(cbind(espectro_PC1, espectro_PC2) ~ forma_PC1 + forma_PC2, 
                     data = pls_df, validation = "LOO")
  
  # Calcular R² e p-valor (como antes)
  predicted <- predict(pls_result, ncomp = 2)
  sse <- sum((pls_df[, c("espectro_PC1", "espectro_PC2")] - predicted)^2)
  sst <- sum(pls_df[, c("espectro_PC1", "espectro_PC2")]^2)
  r_squared <- 1 - (sse/sst)
  
  n_perm <- 999 # Aumentado para resultado mais robusto
  null_r2 <- numeric(n_perm)
  for (i in 1:n_perm) {
    idx_perm <- sample(nrow(pls_df))
    pls_perm <- plsr(cbind(espectro_PC1, espectro_PC2) ~ forma_PC1 + forma_PC2, data = pls_df[idx_perm, ], validation = "none")
    pred_perm <- predict(pls_perm, ncomp = 2)
    sse_perm <- sum((pls_df[, c("espectro_PC1", "espectro_PC2")] - pred_perm)^2)
    null_r2[i] <- 1 - (sse_perm/sst)
  }
  p_value <- mean(null_r2 >= r_squared)
  
  cat(paste("R²:", round(r_squared, 4), "| P-valor:", round(p_value, 4), "\n"))
  
  # --- MUDANÇA IMPORTANTE AQUI ---
  # Retornar uma lista que inclui os scores do PLS para o gráfico
  return(list(
    r_squared = r_squared,
    p_value = p_value,
    scores_df = data.frame(
      especie = common_species,
      forma_pls1 = scores(pls_result)[, 1],    # Scores do eixo X (forma)
      espectro_pls1 = Yscores(pls_result)[, 1] # Scores do eixo Y (espectro)
    )
  ))
}

# =================================================================
# 4. EXECUTAR AS TRÊS ANÁLISES
# =================================================================
cat("\nIniciando análises fPLS...\n")
resultado_geral <- run_phylo_pls_2componentes(contornos_geral, espectros_geral, arvore, "Geral")
resultado_esteril <- run_phylo_pls_2componentes(contorno_esteril, espectro_esteril, arvore, "Estéril")
resultado_fertil <- run_phylo_pls_2componentes(contorno_fertil, espectro_fertil, arvore, "Fértil")
cat("\nAnálises concluídas.\n")

# =================================================================
# 5. CONSTRUÇÃO DA FIGURA
# =================================================================

if (!is.null(resultado_geral) && !is.null(resultado_esteril) && !is.null(resultado_fertil)) {
  
  cat("Construindo a figura conceitual final...\n")
  
  # --- DEFINIÇÃO DA PALETA DE CORES PERSONALIZADA ---
  cores_especies <- c(
    "Microgramma tecta" = "#CC79A7",
    "Microgramma nana" = "#009E73",
    "Microgramma piloselloides" = "#0072B2",
    "Microgramma percussa" = "#F0E442",
    "Microgramma tobagensis" = "#555555",
    "Microgramma latevagans" = "#56B4E9",
    "Microgramma dictyophylla" = "#E69F00",
    "Microgramma reptans" = "#D55E00"
  )
  
  # --- PREPARAÇÃO DE DADOS PARA A FILOGENIA ---
  scores_geral <- resultado_geral$scores_df
  tree_geral <- keep.tip(arvore, tip = scores_geral$especie)
  pls1_forma_vec <- setNames(scores_geral$forma_pls1, scores_geral$especie)
  pls1_espectro_vec <- setNames(scores_geral$espectro_pls1, scores_geral$especie)
  anc_forma <- fastAnc(tree_geral, pls1_forma_vec)
  anc_espectro <- fastAnc(tree_geral, pls1_espectro_vec)
  all_nodes_forma <- c(pls1_forma_vec, anc_forma)
  all_nodes_espectro <- c(pls1_espectro_vec, anc_espectro)
  phylo_segments <- data.frame(
    x_start = all_nodes_forma[tree_geral$edge[, 1]],
    y_start = all_nodes_espectro[tree_geral$edge[, 1]],
    x_end = all_nodes_forma[tree_geral$edge[, 2]],
    y_end = all_nodes_espectro[tree_geral$edge[, 2]]
  )
}
# =================================================================
# 5. CONSTRUÇÃO DA FIGURA NÍVEL NATURE (VERSÃO FINAL COM TAMANHOS MAXIMIZADOS)
# =================================================================

if (!is.null(resultado_geral) && !is.null(resultado_esteril) && !is.null(resultado_fertil)) {
  
  cat("Construindo a figura conceitual final com tamanhos maximizados...\n")
  
  # --- DEFINIÇÃO DA PALETA DE CORES PERSONALIZADA ---
  cores_especies <- c(
    "Microgramma tecta" = "#CC79A7", "Microgramma nana" = "#009E73",
    "Microgramma piloselloides" = "#0072B2", "Microgramma percussa" = "#F0E442",
    "Microgramma tobagensis" = "#555555", "Microgramma latevagans" = "#56B4E9",
    "Microgramma dictyophylla" = "#E69F00", "Microgramma reptans" = "#D55E00"
  )
  
  # --- PREPARAÇÃO DE DADOS PARA A FILOGENIA ---
  scores_geral <- resultado_geral$scores_df
  tree_geral <- keep.tip(arvore, tip = scores_geral$especie)
  pls1_forma_vec <- setNames(scores_geral$forma_pls1, scores_geral$especie)
  pls1_espectro_vec <- setNames(scores_geral$espectro_pls1, scores_geral$especie)
  anc_forma <- fastAnc(tree_geral, pls1_forma_vec)
  anc_espectro <- fastAnc(tree_geral, pls1_espectro_vec)
  all_nodes_forma <- c(pls1_forma_vec, anc_forma)
  all_nodes_espectro <- c(pls1_espectro_vec, anc_espectro)
  phylo_segments <- data.frame(
    x_start = all_nodes_forma[tree_geral$edge[, 1]], y_start = all_nodes_espectro[tree_geral$edge[, 1]],
    x_end = all_nodes_forma[tree_geral$edge[, 2]], y_end = all_nodes_espectro[tree_geral$edge[, 2]]
  )
  
  # --- PAINEL A: O Morfoespaço Evolutivo Geral (TAMANHOS MAXIMIZADOS) ---
  plot_A <- ggplot(scores_geral, aes(x = forma_pls1, y = espectro_pls1)) +
    geom_segment(data = phylo_segments, aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 color = "gray70", alpha = 0.7, linewidth = 1.5) + 
    geom_point(aes(fill = especie), shape = 21, size = 12, color = "black", alpha = 0.9, stroke = 1.5) + # Ponto bem maior
    scale_fill_manual(values = cores_especies, guide = "none") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 8, parse = TRUE, # Anotação maior
             label = sprintf('bold("Combined:")~R^2~"="~%.2f~"|"~italic(P)~"="~%.4f', 
                             resultado_geral$r_squared, resultado_geral$p_value)) +
    labs(title = "", x = "Shape Axis (PLS1)", y = "Spectrum Axis (PLS1)") +
    theme_bw(base_size = 20) + # Base de tamanho aumentada
    theme(plot.title = element_text(face = "bold", size = 24),
          axis.title = element_text(face = "bold", size = 22), # Título do eixo bem maior
          axis.text = element_text(size = 18)) # Texto do eixo bem maior
  
  # --- DADOS PARA O PAINEL B ---
  scores_esteril <- resultado_esteril$scores_df
  scores_fertil <- resultado_fertil$scores_df
  vetores_df <- inner_join(
    scores_esteril %>% rename(forma_esteril = forma_pls1, espectro_esteril = espectro_pls1),
    scores_fertil %>% rename(forma_fertil = forma_pls1, espectro_fertil = espectro_pls1),
    by = "especie"
  )
  
  # --- PAINEL B: O Desacoplamento Visualizado (TAMANHOS MAXIMIZADOS) ---
  plot_B <- ggplot(vetores_df) +
    geom_segment(aes(x = forma_esteril, y = espectro_esteril, xend = forma_fertil, yend = espectro_fertil),
                 color = "black", linewidth = 1.5) + 
    geom_point(aes(x = forma_esteril, y = espectro_esteril, fill = especie), 
               shape = 21, size = 12, color = "black", alpha = 0.9, stroke = 1.5) + # Ponto de partida bem maior
    geom_point(aes(x = forma_fertil, y = espectro_fertil, fill = especie), 
               shape = 21, size = 9, color = "black", alpha = 0.8, stroke = 1.2) + # Ponto de chegada bem maior
    geom_point(aes(x = forma_fertil, y = espectro_fertil), 
               shape = 24, fill = "black", color = "white", size = 4.5, stroke = 0.7) + # Ponta da seta bem maior
    scale_fill_manual(values = cores_especies, name = "") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 8, parse = TRUE, # Anotação maior
             label = sprintf('bold("Sterile:")~R^2~"="~%.2f~"|"~italic(P)~"="~%.4f', 
                             resultado_esteril$r_squared, resultado_esteril$p_value)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 3.5, size = 8, parse = TRUE, # Anotação maior
             label = sprintf('bold("Fertile:")~R^2~"="~%.2f~"|"~italic(P)~"="~%.4f', 
                             resultado_fertil$r_squared, resultado_fertil$p_value)) +
    labs(title = "", subtitle = "",
         x = "Shape Axis (PLS1)", y = "Spectrum Axis (PLS1)") +
    theme_bw(base_size = 20) + # Base de tamanho aumentada
    theme(plot.title = element_text(face = "bold", size = 24),
          axis.title = element_text(face = "bold", size = 22), # Título do eixo bem maior
          axis.text = element_text(size = 18), # Texto do eixo bem maior
          plot.subtitle = element_text(size = 22, face = "italic")) # Subtítulo bem maior
  
  # --- PAINEL C: A Evidência Quantitativa (TAMANHOS MAXIMIZADOS) ---
  df_comparacao <- data.frame(
    Analise = factor(c("Combined", "Sterile", "Fertile"), levels = c("Combined", "Sterile", "Fertile")),
    R_quadrado = c(resultado_geral$r_squared, resultado_esteril$r_squared, resultado_fertil$r_squared),
    Significancia = c("***", "*", "ns")
  )
  cores_analise <- c("Combined" = "#555555", "Sterile" = "#009E73", "Fertile" = "#D55E00")
  plot_C <- ggplot(df_comparacao, aes(x = Analise, y = R_quadrado, fill = Analise)) +
    geom_col(width = 0.7, color = "black", show.legend = FALSE) +
    geom_text(aes(label = Significancia), vjust = -0.5, size = 12) + # Texto de significância bem maior
    scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2), expand = c(0, 0)) +
    scale_fill_manual(values = cores_analise) +
    labs(title = "", x = NULL, y = "Shape and spectrum (R²)") +
    theme_bw(base_size = 20) + # Base de tamanho aumentada
    theme(plot.title = element_text(face = "bold", size = 24),
          axis.title.y = element_text(face = "bold", size = 22), # Título do eixo bem maior
          axis.text.y = element_text(size = 18), # Texto do eixo Y bem maior
          axis.text.x = element_text(size = 20, face = "bold")) # Nomes das barras bem maiores e em negrito
  
  # --- COMBINAR FIGURA FINAL COM LEGENDA NO MEIO ---
  library(cowplot) 
  
  # 1. Extrair a legenda do plot_B (COM TEXTO BEM MAIOR)
  shared_legend <- get_legend(
    plot_B + theme(legend.position = "bottom", 
                   legend.justification = "center",
                   legend.text = element_text(face = "italic", size = 18)) # Legenda bem maior
  )
  
  # 2. Remover a legenda dos plots originais
  plot_A <- plot_A + theme(legend.position = "none")
  plot_B <- plot_B + theme(legend.position = "none")
  
  # 3. Montar o layout de 3 andares
  figura_final_nature <- (plot_A | plot_B) / shared_legend / plot_C +
    plot_layout(heights = c(2.2, 0.3, 1)) # Aumentei um pouco o espaço da legenda
  
  print(figura_final_nature)
  
  ggsave("Imagens_prontas/figura_publicacao_final.png", 
         figura_final_nature, 
         width = 16, height = 12, units = "in", dpi = 300)
  cat("\nFigura final salva como 'figura_publicacao_final.png'\n")
  
} else {
  cat("\nUma ou mais análises falharam. O gráfico não pôde ser gerado.\n")
}

("\n--- FIM DO SCRIPT ---\n")

