# =================================================================
# Script completo para calcular TPM a partir de conteos y longitudes
# Versión final robusta (con prefijos de paquete explícitos)
# =================================================================

# --- Cargar librerías necesarias ---
library(data.table)
library(dplyr)
library(readr)
library(tibble)

# --- Función de utilidad para limpiar IDs de Ensembl (quitar versiones .X) ---
limpiar_ids_ensembl <- function(ids) {
  gsub("\\..*", "", ids)
}
limpiar_nombres <- function(nombres) {
  nombres_acortados <- gsub("^([^.]*\\.[^.]*)\\..*", "\\1", nombres)
  nombres_unicos <- make.unique(nombres_acortados, sep = "_")
  return(nombres_unicos)
}


# --- 1. Leer datos de expresión y metadatos ---
message("1. Leyendo datos de expresión y metadatos...")
leeract_raw <-fread("~/RScripts/genes_expected_count_DESeq2.csv.gz")
ayudaact <- fread("~/RScripts/metadata_monomacro_filtered.csv.gz")

# --- 2. Filtrar metadatos de monocitos específicos ---
message("2. Filtrando metadatos...")
monocyte_meta <- ayudaact %>%
  dplyr::filter(harmonized_donor_id %in% c(
    "S00BS4", "S00C1H", "S00BHQ", "S00FTN", "S001S7", "S0022I",
    "S0018A", "S00DVR", "S007SK", "S00E8W", "S00390", "S006VI",
    "S001MJ", "S00H6O", "S00622", "S00J8C"
  )) %>%
  # ✨ CORRECCIÓN DE ROBUSTEZ: Se usa dplyr::mutate para ser explícito
  dplyr::mutate(across(where(is.character), as.factor))

leeract_raw <- leeract_raw %>%
  rename_with(limpiar_nombres) 

# --- 3. Encontrar muestras comunes entre metadatos y expresión ---
message("3. Identificando muestras comunes...")

common_samples <- intersect(monocyte_meta$EpiRR, colnames(leeract_raw))

# --- 4. Preparando la tabla de conteos ---
message("4. Preparando la tabla de conteos...")
counts <- leeract_raw %>%
  dplyr::select(id_col, all_of(common_samples)) %>%
  dplyr::mutate(id_col = limpiar_ids_ensembl(id_col))

# =================================================================
# ✨ CORRECCIÓN CLAVE: Agregación de conteos y uso de prefijos explícitos
# =================================================================
message("5. Agregando conteos de isoformas para obtener valores por gen...")
counts_aggregated <- counts %>%
  dplyr::group_by(id_col) %>%
  # Se usa dplyr::summarise para evitar conflictos y asegurar compatibilidad
  dplyr::summarise(across(where(is.numeric), sum), .groups = 'drop')
message(paste("   -> Agregación completada. Se obtuvieron", nrow(counts_aggregated), "genes únicos."))
# =================================================================

# --- 6. Leer información genómica y limpiar sus IDs ---
message("6. Leyendo y preparando información genómica (longitudes)...")
nombres <- readr::read_csv("~/RScripts/human_genes_grch37_p13.csv") %>%
  dplyr::mutate(ensembl_gene_id = limpiar_ids_ensembl(ensembl_gene_id))

# --- 7. Diagnóstico: Verificar si hay genes en común ---
message("7. Verificando la coincidencia de genes entre ambos archivos...")
genes_comunes_inicial <- intersect(counts_aggregated$id_col, nombres$ensembl_gene_id)
if (length(genes_comunes_inicial) == 0) {
  stop("❌ ERROR CRÍTICO: No se encontró ninguna coincidencia de genes. Revisa los formatos de los IDs.")
} else {
  message(paste("   -> ✅ Verificación exitosa: Se encontraron", length(genes_comunes_inicial), "genes en común."))
}

# --- 8. Calcular longitudes de genes ---
lengths <- nombres %>%
  dplyr::filter(ensembl_gene_id %in% genes_comunes_inicial) %>%
  dplyr::mutate(lgth = end_position - start_position) %>%
  dplyr::select(ensembl_gene_id, lgth) %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

# --- 9. Crear vector nombrado de longitudes ---
gene_lengths <- lengths$lgth
names(gene_lengths) <- lengths$ensembl_gene_id

# --- 10. Preparar matriz de conteos con rownames ---
message("10. Convirtiendo tabla de conteos a matriz numérica...")
counts_matrix <- as.data.frame(counts_aggregated)
rownames(counts_matrix) <- counts_matrix$id_col
counts_matrix$id_col <- NULL

# ============================================================
# 🧩 11. Función para calcular TPM
# ============================================================
calculate_tpm <- function(counts_mat, gene_len) {
  message("11. Iniciando cálculo de TPM...")
  
  genes_comunes <- intersect(rownames(counts_mat), names(gene_len))
  counts_final <- counts_mat[genes_comunes, , drop = FALSE]
  lengths_final <- gene_len[genes_comunes]
  
  message(paste("    -> Calculando TPM sobre", length(genes_comunes), "genes."))
  
  lengths_kb <- lengths_final / 1000
  rpk <- sweep(counts_final, 1, lengths_kb, "/")
  scaling_factor <- colSums(rpk, na.rm = TRUE) / 1e6
  tpm <- sweep(rpk, 2, scaling_factor, "/")
  
  message("    -> 🎯 TPM calculado correctamente.")
  
  tpm_df <- as.data.frame(tpm) %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    dplyr::as_tibble()
  
  return(list(matrix = as.matrix(tpm), tibble = tpm_df))
}

# --- 12. Ejecutar el cálculo de TPM ---
counts_matrix_filtrada <- counts_matrix[rownames(counts_matrix) %in% names(gene_lengths), ]
resultado_tpm <- calculate_tpm(counts_matrix_filtrada, gene_lengths)

tpm_matrix <- resultado_tpm$matrix
tpm_tibble <- resultado_tpm$tibble


monocyte_meta_filtered <- monocyte_meta %>%
  filter(EpiRR %in% common_samples) %>%
  arrange(match(EpiRR, common_samples))
expr_t <- t(tpm_matrix)

if (!all(rownames(expr_t) == monocyte_meta_filtered$EpiRR)) {
  stop("El orden de las muestras en la matriz de expresión y los metadatos no coincide. ¡Deteniendo!")
} else {
  message("Las muestras en la matriz de expresión y los metadatos están correctamente alineadas.")
}

# --- 4. Cálculo de PCA (UNA SOLA VEZ) ---
message("Realizando cálculo de PCA...")
pca_res <- prcomp(expr_t, center = TRUE, scale. = FALSE)
var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

pca_df_base <- data.frame(
  sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

# --- 5. Función para generar un plot de PCA (SIN LEYENDA en el plot) ---
# Ahora la leyenda se generará aparte para componerla.
generate_pca_plot_core <- function(pca_data_frame, metadata_df, color_column,
                                   title_prefix = "",
                                   point_size = 3, label_points = FALSE,
                                   colors_palette = NULL,
                                   show_legend_for_extraction = FALSE) {
  
  # 1) Hacemos el join (usando dplyr::select con !!sym para seleccionar dinámicamente)
  pca_df_final <- pca_data_frame %>%
    left_join(metadata_df %>% dplyr::select(EpiRR, !!rlang::sym(color_column)),
              by = c("sample" = "EpiRR"))
  
  # 2) Comprobación: existe la columna que queremos usar?
  if (!color_column %in% colnames(pca_df_final)) {
    stop(paste0("La columna '", color_column, "' no existe después del join. ",
                "Revisa el nombre (tiene que coincidir exactamente)."))
  }
  
  # 3) Asignación robusta de la nueva columna color_group (compatible con data.frame y tibble)
  pca_df_final$color_group <- pca_df_final[[color_column]]
  
  # --- resto igual ---
  pc1_var_exp <- round(100 * var_exp[1], 1)
  pc2_var_exp <- round(100 * var_exp[2], 1)
  
  p <- ggplot(pca_df_final, aes(x = PC1, y = PC2, color = color_group)) +
    geom_point(size = point_size) +
    labs(
      title = paste0(title_prefix, "PCA por ", gsub("_", " ", color_column)),
      x = paste0("PC1 (", pc1_var_exp, "%)"),
      y = paste0("PC2 (", pc2_var_exp, "%)"),
      color = gsub("_", " ", color_column)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend_for_extraction) "bottom" else "none",
      legend.justification = "center",
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5, face = "bold")
    )
  
  if (label_points) {
    p <- p + geom_text(aes(label = sample), size = 2.5, vjust = -1, show.legend = FALSE)
  }
  
  if (!is.null(colors_palette)) {
    p <- p + scale_color_manual(values = colors_palette)
  }
  
  return(p)
}

# --- 6. Definir las columnas del metadata por las que quieres iterar ---
message("Identificando columnas adecuadas para colorear los plots...")
excluded_cols_manual <- c("EpiRR", "harmonized_sample_label", "id_col")

columns_to_plot_by <- monocyte_meta_filtered %>%
  dplyr::select(-any_of(excluded_cols_manual)) %>%
  #select(where(~ n_distinct(.) > 1 & n_distinct(.) < (nrow(monocyte_meta_filtered) * 1))) %>%
  colnames()

if (length(columns_to_plot_by) == 0) {
  warning("No se encontraron columnas categóricas adecuadas en los metadatos de monocitos para generar plots.
          Asegúrate de tener columnas con al menos dos categorías y no demasiadas.")
} else {
  message(paste("Se generarán plots para las siguientes columnas:", paste(columns_to_plot_by, collapse = ", ")))
}

# --- 7. Generar y guardar plots con leyenda separada debajo en un único PDF (tamaño A4) ---
composed_plots_list <- list() # Lista para almacenar los plots compuestos
output_dir <- "PCA_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
output_pdf_path <- file.path(output_dir, "MacroTPM_Plots_A4.pdf")

color_palettes_base <- c("Set1", "Dark2", "Paired", "Accent", "Spectral", "RdYlGn", "BrBG")

if (length(columns_to_plot_by) > 0) {
  a4_width_inches <- 8.27
  a4_height_inches <- 11.69

  message(paste("Guardando plots en:", output_pdf_path, "con tamaño A4."))

  pdf(output_pdf_path, width = a4_width_inches, height = a4_height_inches)

  for (col_name in columns_to_plot_by) {
    message(paste("Procesando columna para plot:", col_name))
    monocyte_meta_filtered[[col_name]] <- as.factor(monocyte_meta_filtered[[col_name]])
    num_levels <- nlevels(monocyte_meta_filtered[[col_name]])

    # Lógica de selección de paleta de colores para mayor contraste
    if (num_levels <= 8) { # Para un número pequeño de categorías, usamos Dark2 o Set1
      current_colors <- RColorBrewer::brewer.pal(n = max(3, num_levels), name = "Dark2")[1:num_levels]
    } else if (num_levels <= 12) { # Para un número intermedio, Paired suele funcionar bien
      current_colors <- RColorBrewer::brewer.pal(n = max(3, num_levels), name = "Paired")[1:num_levels]
    } else { # Para muchos niveles, usamos una rampa de color más amplia y contrastada
      # Puedes ajustar las paletas base aquí para experimentar con diferentes rangos
      current_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_levels)
    }

    # 1. Generar el plot que usaremos para extraer la leyenda
    plot_for_legend_extraction <- generate_pca_plot_core(
      pca_data_frame = pca_df_base,
      metadata_df = monocyte_meta_filtered,
      color_column = col_name,
      title_prefix = "",
      point_size = 2.5,
      label_points = FALSE,
      colors_palette = current_colors,
      show_legend_for_extraction = TRUE # Genera la leyenda en la parte inferior para que get_legend la encuentre
    )

    # 2. Extraer la leyenda
    plot_legend <- get_legend(plot_for_legend_extraction)

    # 3. Generar el plot final SIN leyenda (ya que la vamos a componer aparte)
    main_plot <- generate_pca_plot_core(
      pca_data_frame = pca_df_base,
      metadata_df = monocyte_meta_filtered,
      color_column = col_name,
      title_prefix = "",
      point_size = 2.5,
      label_points = FALSE,
      colors_palette = current_colors,
      show_legend_for_extraction = FALSE # Sin leyenda en el plot principal
    )

    # 4. Componer el plot principal y la leyenda usando plot_grid
    # main_plot toma la mayor parte del espacio, legend_plot una fracción menor
    # rel_heights controla la proporción (ej. 0.8 para el plot, 0.2 para la leyenda)
    composed_plot <- plot_grid(
      main_plot,
      plot_legend,
      ncol = 1,
      rel_heights = c(0.8, 0.2), # Ajusta estos valores para controlar el espacio de la leyenda
      align = "v" # Alinea verticalmente
    )

    # Imprimir el plot compuesto en el PDF
    print(composed_plot)
  }

  dev.off() # Cierra el dispositivo PDF
  message("Plots con leyenda individual y tamaño controlado guardados exitosamente en PDF.")

} else {
  message("No se generaron plots debido a la falta de columnas válidas.")
}

message("Proceso completo.")
TPMSUM<- as.data.frame(t(colSums(tpm_matrix)))
EXPSUM<-format(as.data.frame(t(colSums(counts[ , -1], na.rm = TRUE))),scientific=TRUE)

pc1_values <- pca_res$x[, 1]
names(pc1_values) <- rownames(pca_res$x)

# Identificar las muestras con PC1 > 0 y PC1 < 0
samples_pc1_pos <- names(pc1_values[pc1_values > 0])
samples_pc1_neg <- names(pc1_values[pc1_values < 0])

# Filtrar EXPSUM o TPMSUM según estas muestras
# (Nota: EXPSUM tiene las columnas en el mismo orden que las muestras)
EXPSUM_PC1_pos <- EXPSUM %>%
  dplyr::select(all_of(samples_pc1_pos)) %>%
  dplyr::mutate(grupo_PC1 = "positivo")

EXPSUM_PC1_neg <- EXPSUM %>%
  dplyr::select(all_of(samples_pc1_neg)) %>%
  dplyr::mutate(grupo_PC1 = "negativo")

