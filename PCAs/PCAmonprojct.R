library(readr)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(data.table)

leeract_raw<-fread("~/RScripts/genes_tpm_monomacro_processed.csv.gz")
ayudaact<-fread("~/RScripts/metadata_monomacro_filtered.csv.gz")
monocyte<- ayudaact %>%
  filter(harmonized_sample_label == "Monocyte")
macrophage<- ayudaact %>%
  filter(harmonized_sample_label == "Macrophage")

monocyte_proyect <- split(monocyte,monocyte$project)

colsmon <- c("id_col",intersect(monocyte$EpiRR, colnames(leeract_raw)))

act2 <- leeract_raw[, ..colsmon]

expr_matrix <- act2


rownames(expr_matrix) <- expr_matrix$id_col

expr_matrix <- expr_matrix[ , -1]
expr_t <- t(expr_matrix)

pca <- prcomp(expr_t, center = TRUE, scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)



# Data frame con los resultados de la PCA
pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

# Marcar si la muestra está en tu lista/meta
pca_df$color <- ifelse(pca_df$sample %in% monocyte_proyect$BLUEPRINT$EpiRR, "blue",
                       ifelse(pca_df$sample %in% monocyte_proyect$CEEHRC$EpiRR, "red",
                              ifelse(pca_df$sample %in% monocyte_proyect$DEEP$EpiRR, "green",
                                     ifelse(pca_df$sample %in% monocyte_proyect$ENCODE$EpiRR, "orange",
                                            ifelse(pca_df$sample %in% monocyte_proyect$`NIH Roadmap Epigenomics`$EpiRR, "purple", "grey")))))

# Plot
plot(pca_df$PC1, pca_df$PC2,
     col = pca_df$color,
     pch = 19,
     xlab = paste0("PC1 (", round(100*var_exp[1],1), "%)"),
     ylab = paste0("PC2 (", round(100*var_exp[2],1), "%)"),
     main = "PCA de los monocitos por proyecto")

#legend("topright", legend = c("Monocyte", "Macrophage"),
     #  col = c("blue", "red"), pch = 19)

#text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample,     pos = 3, cex = 0.7)

