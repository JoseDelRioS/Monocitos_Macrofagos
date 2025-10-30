library(readr)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(cowplot) # Fundamental para get_legend y plot_grid
library(RColorBrewer)
library(grid) # Para grid.text si queremos añadir títulos personalizados

# --- 1. Carga de datos ---
message("Cargando datos...")
leeract_raw <- fread("~/RScripts/genes_expected_count_DESeq2.csv.gz")
ayudaact <- fread("~/RScripts/metadata_monomacro_filtered.csv.gz")


# --- 2. Preparación inicial de metadatos ---
monocyte_meta <- ayudaact %>%
  filter(harmonized_donor_id %in% c(
    "S00BS4", "S00C1H", "S00BHQ", "S00FTN", "S001S7", "S0022I",
    "S0018A", "S00DVR", "S007SK", "S00E8W", "S00390", "S006VI",
    "S001MJ", "S00H6O", "S00622", "S00J8C"
  )) %>%
  mutate(across(where(is.character), as.factor))
# --- 3. Preparación de la matriz de expresión ---
message("Preparando matriz de expresión...")

limpiar_nombres <- function(nombres) {
  nombres_acortados <- gsub("^([^.]*\\.[^.]*)\\..*", "\\1", nombres)
  nombres_unicos <- make.unique(nombres_acortados, sep = "_")
  return(nombres_unicos)
}

leeract_raw <- leeract_raw %>%
  rename_with(limpiar_nombres) 

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
                                   show_legend_for_extraction = FALSE) { # Para extraer la leyenda o no
  
  pca_df_final <- pca_data_frame %>%
    left_join(metadata_df %>% select(EpiRR, !!sym(color_column)),
              by = c("sample" = "EpiRR")) %>%
    rename(color_group = !!sym(color_column))
  
  pc1_var_exp <- round(100*var_exp[1], 1)
  pc2_var_exp <- round(100*var_exp[2], 1)
  
  p <- ggplot(pca_df_final, aes(x = PC1, y = PC2, color = color_group)) +
    geom_point(size = point_size) +
    labs(
      title = paste0(title_prefix, "PCA por ", gsub("_", " ", color_column)),
      x = paste0("PC1 (", pc1_var_exp, "%)"),
      y = paste0("PC2 (", pc2_var_exp, "%)"),
      color = gsub("_", " ", color_column) # Leyenda del color
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend_for_extraction) "bottom" else "none", # CAMBIO CLAVE
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
  select(-any_of(excluded_cols_manual)) %>%
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
output_pdf_path <- file.path(output_dir, "MonocitosRAw_Plots_A4.pdf")

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