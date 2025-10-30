library(readr)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(grid)
library(forcats)

# --- 1. Carga de datos ---
message("Cargando datos...")
leeract_raw <- fread("~/RScripts/genes_tpm_monomacro_processed.csv.gz")
ayudaact <- fread("~/RScripts/metadata_monomacro_filtered.csv.gz")

# --- 2. Preparación inicial de metadatos ---
monocyte_meta <- ayudaact %>%
filter(harmonized_sample_label == "Monocyte") %>%
 mutate(across(where(is.character), as.factor))
#Doble columna

# monocyte_meta <- ayudaact %>%
#   filter(harmonized_sample_label == "Macrophage") %>%
#   as.data.frame() %>%                          # <- lo pasamos a data.frame
#   mutate(across(everything(), as.character)) %>% 
#   {
#     df <- .
#     
#     # excluir la primera columna
#     cols <- names(df)[-1]
#     
#     # quedarnos solo con columnas NO constantes
#     cols_var <- cols[sapply(df[cols], function(x) length(unique(x)) > 1)]
#     
#     # todas las combinaciones posibles de 2 en 2
#     combos <- combn(cols_var, 2, simplify = FALSE)
#     
#     # generar nuevas columnas
#     new_cols <- lapply(combos, function(cols) {
#       paste(df[[cols[1]]], df[[cols[2]]])
#     })
#     
#     # nombres de las nuevas columnas
#     names(new_cols) <- sapply(combos, function(cols) paste(cols, collapse = "_"))
#     
#     cbind(df, as.data.frame(new_cols))
#   }%>%
#   mutate(across(where(is.character), as.factor))

# --- 3. Preparación de la matriz de expresión ---
message("Preparando matriz de expresión...")
common_samples <- intersect(monocyte_meta$EpiRR, colnames(leeract_raw))

expr_matrix_raw <- leeract_raw %>%
  select(id_col, all_of(common_samples))

monocyte_meta_filtered <- monocyte_meta %>%
  filter(EpiRR %in% common_samples) %>%
  arrange(match(EpiRR, common_samples))

expr_matrix <- as.data.frame(expr_matrix_raw)
rownames(expr_matrix) <- expr_matrix$id_col
expr_matrix <- expr_matrix[, -1]
expr_matrix <- as.matrix(expr_matrix)

expr_t <- t(expr_matrix)

if (!all(rownames(expr_t) == monocyte_meta_filtered$EpiRR)) {
  stop("El orden de las muestras en la matriz de expresión y los metadatos no coincide. ¡Deteniendo!")
} else {
  message("Las muestras en la matriz de expresión y los metadatos están correctamente alineadas.")
}

# --- 4. PCA ---
message("Realizando cálculo de PCA...")
pca_res <- prcomp(expr_t, center = TRUE, scale. = FALSE)
var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

pca_df_base <- data.frame(
  sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

pca_df_con_meta <- pca_df_base %>%
  rename(EpiRR = sample) %>% 
  left_join(monocyte_meta_filtered, by = "EpiRR")

# --- 5. Añadimos la columna 'categoria' ---
pca_df_final <- pca_df_con_meta %>%
  mutate(categoria = if_else(PC1 < -1e+05, "A", "B"))

# --- 6. Preparamos los datos ---
datos_ml <- pca_df_final %>%
  select(-PC1, -PC2,-EpiRR_ordering,-EpiRR) %>%
  mutate(
    across(where(is.character), as.factor),
    across(where(is.factor), ~ fct_explicit_na(., na_level = "Missing")),
    categoria = as.factor(categoria)
  )

# --- 7. Inicializamos H2O ---
library(h2o)
h2o.init(nthreads = -1, max_mem_size = "8G")

# Convertimos a formato H2O
datos_h2o <- as.h2o(datos_ml)

# Definimos variables predictoras y target
y <- "categoria"
x <- setdiff(names(datos_ml), y)

# --- 8. Dividimos en train/test ---
splits <- h2o.splitFrame(datos_h2o, ratios = 0.8, seed = 123)
train_h2o <- splits[[1]]
test_h2o  <- splits[[2]]

# --- 9. Entrenamos red neuronal ---
dl_model <- h2o.deeplearning(
  x = x,
  y = y,
  training_frame = train_h2o,
  validation_frame = test_h2o,
  activation = "RectifierWithDropout",
  hidden = c(64, 32),   # capas ocultas
  epochs = 500,
  l1 = 1e-5,
  l2 = 1e-5,
  seed = 123,
  balance_classes = TRUE
)

# --- 10. Evaluación ---
sh <- h2o.scoreHistory(dl_model)
perf <- h2o.performance(dl_model, test_h2o)
print(perf)



# --- 11. Importancia de variables ---
var_imp <- h2o.varimp(dl_model)
print(var_imp)
h2o.varimp_plot(dl_model)
