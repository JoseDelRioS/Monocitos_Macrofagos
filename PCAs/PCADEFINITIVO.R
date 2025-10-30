library(readr)
library(GenomicRanges) # Aunque no se usa en este código, lo mantengo por si lo necesitas más adelante.
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2) # Para gráficos más bonitos y flexibles
library(cowplot) # Para organizar múltiples plots si los generas en un loop

# --- 1. Carga de datos (mantengo tus rutas) ---
# Usar data.table para cargar es eficiente.
message("Cargando datos...")
leeract_raw <- fread("~/RScripts/genes_tpm_monomacro_processed.csv.gz")
ayudaact <- fread("~/RScripts/metadata_monomacro_filtered.csv.gz")

# --- 2. Preparación inicial de metadatos ---
# Puedes combinar los filtros iniciales si es que los necesitas por separado.
# Por ahora, los mantengo como los tienes, pero considera si realmente necesitas
# objetos 'monocyte' y 'macrophage' separados si vas a iterar sobre 'ayudaact'.

# Si vas a iterar sobre todas las columnas de 'ayudaact', tal vez no necesites
# los filtros iniciales de 'Monocyte' y 'Macrophage' de forma rígida al principio.
# Podrías querer que el loop se haga para 'ayudaact' completo o para un subconjunto de él.

# Para el ejemplo, asumiré que quieres iterar sobre las columnas de 'monocyte'
# y aplicar la PCA a 'monocyte' específicamente.
monocyte_meta <- ayudaact %>%
  filter(harmonized_sample_label == "Monocyte")

macrophage_meta <- ayudaact %>%
  filter(harmonized_sample_label == "Macrophage")

# --- 3. Preparación de la matriz de expresión ---
message("Preparando matriz de expresión...")
# Solo necesitamos las columnas de expresión que están en el metadata de monocyte
# y la columna de ID.
# Usamos intersect para asegurarnos de que solo seleccionamos columnas existentes.
common_samples <- intersect(monocyte_meta$EpiRR, colnames(leeract_raw))
# Crear colsmon con 'id_col' y las muestras comunes
colsmon <- c("id_col", common_samples)
expr_matrix_raw <- leeract_raw[, ..colsmon]

# Establecer id_col como rownames y remover la columna id_col
expr_matrix <- as.data.frame(expr_matrix_raw) # Convertir a data.frame para rownames más fácilmente
rownames(expr_matrix) <- expr_matrix$id_col
expr_matrix <- expr_matrix[ , -1] # Eliminar la columna id_col
expr_matrix <- as.matrix(expr_matrix) # Asegurarse de que sea una matriz numérica

# Transponer la matriz para PCA (muestras como filas, genes como columnas)
expr_t <- t(expr_matrix)

# Asegurarse de que las muestras en expr_t y monocyte_meta estén en el mismo orden
# Esto es CRÍTICO para la asignación correcta de colores.
monocyte_meta_ordered <- monocyte_meta %>%
  filter(EpiRR %in% rownames(expr_t)) %>%
  arrange(match(EpiRR, rownames(expr_t)))

# Verificar si el orden es el mismo (debería ser TRUE)
all(rownames(expr_t) == monocyte_meta_ordered$EpiRR)

# --- 4. Función para realizar PCA y generar el plot ---
# Encapsulamos la lógica en una función para reutilizarla fácilmente.
generate_pca_plot <- function(metadata_df, expr_transposed, color_column, title_prefix = "") {
  message(paste("Generando PCA para:", color_column))
  
  # 1. Realizar PCA
  pca_res <- prcomp(expr_transposed, center = TRUE, scale. = FALSE)
  var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
  
  # 2. Crear data frame con los resultados de la PCA
  pca_df <- data.frame(
    sample = rownames(pca_res$x),
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2]
  )
  
  # 3. Añadir la columna de color del metadata
  # Asegúrate de que las muestras coincidan entre pca_df y metadata_df
  # Usa un left_join para unir los datos por la columna 'sample'/'EpiRR'
  pca_df_final <- pca_df %>%
    left_join(metadata_df %>% select(EpiRR, !!sym(color_column)),
              by = c("sample" = "EpiRR")) %>%
    # Renombrar la columna a 'color_group' para ser más genérico en el plot
    rename(color_group = !!sym(color_column))
  n_groups <- length(unique(pca_df_final$color_group))
  # 4. Generar el plot con ggplot2
  p <- ggplot(pca_df_final, aes(x = PC1, y = PC2, color = color_group)) +
    geom_point(size = 3) +
    labs(
      title = paste0(title_prefix, "PCA por ", gsub("_", " ", color_column)), # Mejora el título
      x = paste0("PC1 (", round(100*var_exp[1],1), "%)"),
      y = paste0("PC2 (", round(100*var_exp[2],1), "%)"),
      color = gsub("_", " ", color_column) # Leyenda del color
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.6 + 0.05 * n_groups, "cm") 
    )
  
  return(p)
}

# --- 5. Definir las columnas del metadata por las que quieres iterar ---
# Excluye columnas que no tienen sentido para colorear (ej. IDs, columnas numéricas continuas que no son factores)
# Ojo: Si tienes columnas con muchos valores únicos (ej. IDs), esto generará muchos colores y será ilegible.
# Mejor elige columnas categóricas o discretas.
# Puedes inspeccionar monocyte_meta para decidir qué columnas usar:
# str(monocyte_meta)
# unique(monocyte_meta$harmonized_donor_sex)
# unique(monocyte_meta$harmonized_ethnicity) # Si hay muchas, considera agruparlas.

# Columnas de monocyte_meta que son de interés para la coloración
# Asegúrate de que estas columnas sean factores o se puedan tratar como tales.
# Por ejemplo, si una columna numérica tiene pocos valores únicos, ggplot la manejará bien.
# Obtener todos los nombres de las columnas
all_names <- names(monocyte_meta)

# Crear un nuevo vector excluyendo el primer nombre
columns_to_plot_by <- all_names[-1]

# Puedes imprimirlo para verificar


#columns_to_plot_by <- c(
#  "project",
#  "harmonized_biomaterial_type" ,
#  "harmonized_sample_ontology_term_high_order_fig1_color" ,
#  "automated_harmonized_sample_ontology"  
#)

# Filtra solo las columnas que realmente existen en monocyte_meta
columns_to_plot_by <- intersect(columns_to_plot_by, colnames(monocyte_meta))

# Filtra solo las columnas que realmente existen en monocyte_meta
columns_to_plot_by <- intersect(columns_to_plot_by, colnames(monocyte_meta))

# --- 6. Iterar y generar plots individualmente ---
# NO ALMACENAMOS EN UNA LISTA SI SOLO VAMOS A MOSTRARLAS INDIVIDUALMENTE,
# SIMPLEMENTE LAS IMPRIMIMOS UNA POR UNA.
if (length(columns_to_plot_by) > 0) {
  for (col_name in columns_to_plot_by) {
    # Asegúrate de que la columna se trate como factor
    monocyte_meta_ordered[[col_name]] <- as.factor(monocyte_meta_ordered[[col_name]])
    
    current_plot <- generate_pca_plot(
      metadata_df = monocyte_meta_ordered,
      expr_transposed = expr_t,
      color_column = col_name,
      title_prefix = "Monocyte "
    )
    
    print(current_plot) # Imprime cada plot individualmente
    # Aquí puedes añadir código para guardar cada plot si lo deseas
    # ggsave(filename = paste0("pca_monocyte_", col_name, ".png"), plot = current_plot, width = 8, height = 7)
    # ggsave(filename = paste0("pca_monocyte_", col_name, ".pdf"), plot = current_plot, width = 8, height = 7)
  }
} else {
  message("No se encontraron columnas válidas para generar plots. Por favor, verifica 'columns_to_plot_by'.")
}

message("Proceso de generación de plots individuales terminado.")