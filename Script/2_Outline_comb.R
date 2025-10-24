# =============================================================================
# OUTLINE_comb
# =============================================================================
# Article: ........................... 
# Authors: Niksoney A. Mendonça, Juliana Aljahara, D. Victor S. Silva and Thaís E. Almeida
# Script author: Niksoney Azevedo Mendonça
# E-mail: niksoneyazevedo2017@gmail.com
# =============================================================================
# Instalar e Carregar Pacotes
# =============================================================================
# Instalar todos os pacotes de uma vez
install.packages(c("Momocs","readxl","tidyverse","ggplot2","tidyr","ggExtra","cowplot","dplyr","ggrepel"))

# Carregar todos os pacotes
library(Momocs)    # Análise de contornos e formas (morfometria geométrica)
library(readxl)    # Ler arquivos Excel
library(tidyverse) # Coleção de pacotes integrados para análise de dados
library(ggplot2)   # Criação de gráficos elegantes
library(tidyr)     # Organizar dados em formato limpo
library(ggExtra)   # Adicionar gráficos marginais ao ggplot2
library(cowplot)   # Combinar e alinhar gráficos ggplot2
library(dplyr)     # Manipulação de dados com sintaxe intuitiva
library(ggrepel)   # Melhorar a legibilidade de rótulos em gráficos

# =============================================================================
# Definindo 8 cores manualmente que vou usar
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7", "#555555")

#-----------------------------------------------------------------------------------------------------------------------
# Processamento inicial dos dados
#-----------------------------------------------------------------------------------------------------------------------

# Definir o diretório de trabalho para onde estão os arquivos
# Colocar o nome da pasta da espécie que to trabalhando
# Ler os arquivos de texto na pasta

# Configurar diretório de trabalho
setwd("C:/Users/nikso/OneDrive/Vida acadêmica e pessoal - Niksoney/2_Biologia Vegetal_UFPE (Mestrado_2023-2024)/_NIKSONEY AZEVEDO_mestrado/2 - MORFOMETRIA GEOMÉTRICA - (MG)/_Oultine")

# Listar arquivos de texto para ambos os tipos de folhas
lf_esteril <- list.files("Coordenadas_fourier_estéril/", pattern = "\\.txt$", full.names = TRUE)
lf_fertil <- list.files("Coordenadas_fourier_fértil/", pattern = "\\.txt$", full.names = TRUE)

# Ler os dados de ambos os tipos
spp_names_esteril <- read_xlsx("Coordenadas_fourier_estéril/tabelaoutline_esteril.xlsx")
spp_names_fertil <- read_xlsx("Coordenadas_fourier_fértil/tabelaoutline_fertil.xlsx")

# Adicionar uma coluna indicando o tipo de folha
spp_names_esteril$tipo <- "estéril"
spp_names_fertil$tipo <- "fértil"

# Combinar os dados
spp_names_combinado <- rbind(spp_names_esteril, spp_names_fertil)

# Combinar a lista de arquivos
lf_combinado <- c(lf_esteril, lf_fertil)

# Exibir os nomes únicos das espécies na coluna 'spp' do arquivo Excel
unique(spp_names_combinado$especie)

# Importar as coordenadas dos arquivos de texto para um objeto de dados
lf_coord <- import_txt(lf_combinado)

# Verificar se o número de arquivos e linhas na planilha são iguais
stopifnot(length(lf_combinado) == nrow(spp_names_combinado))

# Renomear arquivos com base no identificador
names(lf_combinado) <- spp_names_combinado$individuo

# Atribuir os nomes das espécies aos objetos
names(lf_coord) <- spp_names_combinado$especie

# Criar data frame com fatores
lf_fac3 <- data.frame(
  Type = spp_names_combinado$especie,
  Heb = spp_names_combinado$herbário
)

# Criar objeto do Momocs
lf_out3 <- Out(lf_coord, fac = lf_fac3)
lf_out3$fac$Type <- as.factor(lf_out3$fac$Type)
#-----------------------------------------------------------------------------------------------------------------------
# Formas dos grupos
#-----------------------------------------------------------------------------------------------------------------------

# Primeiro abre o dispositivo gráfico PNG
png("Imagens_prontas/formas_combinadas.png", width = 10, height = 8, units = "in", res = 300)

# Depois gera o gráfico
formas_gerais <- panel(lf_out, 
                       fac = "Type", 
                       names = "",
                       col = my_colors)

# Adiciona a legenda
legend("topright", 
       legend = levels(lf_out$Type),
       col = my_colors[1:length(unique(lf_out$Type))],
       pch = 16,
       pt.cex = 2.5,
       cex = 0.8,
       bty = "n",
       inset = c(0, 0),
       xpd = TRUE,
       y.intersp = 2)  # Aumenta o espaçamento vertical entre os itens

# Finalmente fecha o dispositivo gráfico
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# Visualização de contornos
#-----------------------------------------------------------------------------------------------------------------------

stack(coo_center(lf_out),           #Centraliza os contornos antes de plotá-los
      palette = col_summer,         #Define a paleta de cores para o gráfico
      fac = lf_out$Type,            #Fator que define como os contornos são agrupados ou coloridos
      title = "Stack of outlines from all individuals",
      subtitle = paste("Relação entre Área e Comprimento dos Contornos das Folhas"))  #Adiciona um subtítulo com informações extras

#-----------------------------------------------------------------------------------------------------------------------
# PCA e LDA com dados sem processamento_EFA
#-----------------------------------------------------------------------------------------------------------------------

# Realizando a análise de Fourier elíptica (EFA) nas formas armazenadas em lf_out
lf_fou3 <- efourier(x = lf_out3,    #Dados das coordenadas das formas (contornos) 
                   nb.h = 9,     #Quantos harmonicos o comando anterior falou?
                   norm = TRUE)  #Evitar normalizar as formas (norm = FALSE)

#----------------------PCA

# Realizar a PCA
lf_pca3 <- PCA(x = lf_fou3,        #Objeto com os coeficientes de Fourier
              fac = lf_fou3$fac)  #Variável categórica associada (fac)

#-----------------------------------------------------------------------------------------------------------------------
# Função PCA com centroides destacados e marginais (sem elipses)
#-----------------------------------------------------------------------------------------------------------------------
criar_plot_pca <- function(pca_obj, out_obj) {
  pca_data <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    especie = out_obj$fac$Type
  )
  
  # % de variância
  variance_explained <- round((pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100, 1)
  
  # Centroides
  centroids <- pca_data %>%
    group_by(especie) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2))
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = especie)) +
    geom_point(size = 3, alpha = 0.7) +
    
    # CENTRÓIDES → maiores, com borda preta
    geom_point(data = centroids, aes(x = PC1, y = PC2, fill = especie),
               color = "black", shape = 21, size = 6, stroke = 1.2) +
    
    # RÓTULOS nos centroides
    geom_text_repel(data = centroids, aes(x = PC1, y = PC2, label = especie),
                    color = "black", size = 4, fontface = "bold", 
                    box.padding = 0.3, point.padding = 0.3, show.legend = FALSE) +
    
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    labs(
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      axis.title.x = element_text(size = 14, face = "bold"),  # ↑ tamanho maior e negrito
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  # Marginais
  p_marginal <- ggMarginal(
    p, 
    type = "density", 
    groupColour = TRUE,
    groupFill = TRUE
  )
  return(p_marginal)
}

#-----------------------------------------------------------------------------------------------------------------------
# Gerar gráficos
#-----------------------------------------------------------------------------------------------------------------------
p3 <- criar_plot_pca(lf_pca3, lf_out3)

# Salvar
ggsave("Imagens_prontas/pca_combinada_ggplot.png", figura_final,
       width = 12, height = 15, dpi = 300)

# Exibir
print(p3)

#----------------------LDA

# Definir o nome do arquivo e a resolução (300 DPI)
png("Imagens_prontas/lda_combinado.png", width = 10, height = 8, units = "in", res = 300)

# Dados já existentes (supondo que lf_fou já foi criado antes da PCA)
coefs <- lf_fou$coe           # Extrai coeficientes
fac_data <- lf_fou$fac        # Extrai fatores (Type e Heb)

# Calcular MÉDIAS por Heb (indivíduo)
lf_fou_avg <- data.frame(coefs, fac_data) %>%  
  group_by(Heb, Type) %>%  
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Remover as colunas com valores 0 e 1
lf_fou_avg <- lf_fou_avg[, -c(3, 12, 21)]

#3, 9, 15 - 6 harmonicos
#3, 10, 17 - 7 harmonicos
#3, 11, 19 - 8 harmonicos
#3, 12, 21 - 9 harmonicos
#3, 13, 23 - 10 harmonicos

# Converter para formato do Momocs e rodar LDA
lda_data1 <- LDA(
  x = as.matrix(lf_fou_avg[, -c(1:2)]),  # Exclui colunas Heb e Type
  fac = lf_fou_avg$Type)                  # Usa Type como grupo

# Plotar LDA
plot_LDA(lda_data1, palette = pal_manual(my_colors, transp = 0))

# Mostrar resultados (incluindo acurácia)
print(lda_data1)  # Mostra estatísticas gerais

# Acurácia específica:
cat("\nAcurácia da LDA:", mean(lda_data1$CV.correct) * 100, "%\n")

# Matriz de confusão
cat("\nMatriz de Confusão (Validação Cruzada):\n")
print(lda_data1$CV.tab)

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

############################
# matriz de confusão tipo 1#
############################

# Definir o nome do arquivo e a resolução (300 DPI)
png("Imagens_prontas/cv_combinado.png", width = 10, height = 8, units = "in", res = 300)

# Matriz de confusão
plot_CV(lda_data1,
        freq = TRUE,
        rm0 = TRUE,
        fill = TRUE,
        axis.size = 15,
        axis.x.angle = 50,
        cell.size = 6)

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

############################
# matriz de confusão tipo 2#
############################

# Definir o nome do arquivo e a resolução (300 DPI)
png("Imagens_prontas/cv2_combinado.png", width = 10, height = 8, units = "in", res = 300)

plot_CV2(lda_data1,
         freq = FALSE,
         rm0 = TRUE,
         fill = TRUE,
         axis.size = 15,
         axis.x.angle = 50,
         cell.size = 6)

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# Plotar as médias dos formatos das espécies após a LDA
#-----------------------------------------------------------------------------------------------------------------------

# Definir o nome do arquivo e a resolução (300 DPI)
png("Imagens_prontas/formas_médiascv_combinado.png", width = 10, height = 8, units = "in", res = 300)

# Definir as cores
my_colors_shape <- c("#0049ff","#FF7F50")

# Calcular as médias das formas
lf_mean_shape <- MSHAPES(lf_fou, ~Type)

# Plotar usando plot_MSHAPES (sem o argumento 'main')
plot_MSHAPES(lf_mean_shape, palette = pal_manual(my_colors_shape, transp = 0))

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# Análise multivariada da variância - MANOVA
#-----------------------------------------------------------------------------------------------------------------------

# Rodar a MANOVA para o modelo
MANOVA(lf_pca, ~Type)

# Rodar a MANOVA par a par por especie
resultado_manova <- MANOVA_PW(lf_pca, ~Type)  

# Salvar a tabela de significância ($stars.tab) em um arquivo CSV
write.csv(resultado_manova$stars.tab, file = "Imagens_prontas/manova_stars_combinado.csv", row.names = TRUE)

dados_manova <- read.csv("Imagens_prontas/manova_stars_combinado.csv", sep = ",", row.names = 1)

# Transformação dos dados incluindo o nível marginal
dados_long <- dados_manova %>%
  mutate(Observação = row.names(dados_manova)) %>%
  pivot_longer(
    cols = -Observação,
    names_to = "Variável",
    values_to = "Significância"
  ) %>%
  mutate(
    Significância_desc = case_when(
      Significância == "*" ~ "italic('p') < 0.05",
      Significância == "**" ~ "italic('p') < 0.01",
      Significância == "***" ~ "italic('p') < 0.001",
      Significância == "." ~ "italic('p') < 0.1",  # Novo nível adicionado
      TRUE ~ NA_character_
    )
  )

# Definir o nome do arquivo e a resolução (300 DPI)
png("Imagens_prontas/MANOVA_par_a_par_combinado.png", width = 10, height = 8, units = "in", res = 300)

# Heatmap atualizado com a nova categoria
ggplot(dados_long, aes(x = Variável, y = Observação, fill = Significância_desc)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    name = "Significance Level",
    values = c(
      "italic('p') < 0.1" = "#ffffa1", 
      "italic('p') < 0.05" = "gold",  
      "italic('p') < 0.01" = "darkorange",  
      "italic('p') < 0.001" = "chocolate"   
    ),
    na.value = "grey90",
    labels = c(
      "italic('p') < 0.1" = expression(italic('p') < 0.1),
      "italic('p') < 0.05" = expression(italic('p') < 0.05),
      "italic('p') < 0.01" = expression(italic('p') < 0.01),
      "italic('p') < 0.001" = expression(italic('p') < 0.001)
    ),
    drop = TRUE,
    na.translate = FALSE
  ) +
  theme_minimal() +
  labs(x = "Variable", y = "Observation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# Plotas as formas médias por espécie
#-----------------------------------------------------------------------------------------------------------------------

# 1. Calcular as formas médias
formas_medias <- MSHAPES(lf_fou, ~Type)

# 2. Verificar a estrutura do objeto
print(str(formas_medias))  # Isso mostra como o objeto está organizado

# 3. Criar objeto Out corretamente
# O objeto MSHAPES retorna as formas em $shp
formas_medias_out <- Out(formas_medias$shp)

# 4. Determinar o número de grupos
num_grupos <- length(formas_medias$shp)  # Alternativa: nlevels(lf_fou$fac$Type)

# Definir o nome do arquivo e a resolução (300 DPI)
png("Imagens_prontas/formas_ind_combinada.png", width = 10, height = 8, units = "in", res = 300)

# 5. Plotar com cores personalizadas
panel(formas_medias_out, 
      border = my_colors[1:num_grupos],  # Contour colors
      col = NA,                          # No fill
      lwd = 8,                           # Line width
      names = "")                     # Correct way to suppress names

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()