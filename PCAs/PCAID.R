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
leeract_raw <- fread("~/RScripts/genes_tpm_monomacro_processed.csv.gz")
ayudaact <- fread("~/RScripts/metadata_monomacro_filtered.csv.gz")


# --- 2. Preparación inicial de metadatos ---
monocyte_meta <- ayudaact %>%
  filter(harmonized_sample_label == "Macrophage") %>%
  mutate(across(where(is.character), as.factor))

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

# --- 4. Cálculo de PCA (UNA SOLA VEZ) ---
message("Realizando cálculo de PCA...")
pca_res <- prcomp(expr_t, center = TRUE, scale. = FALSE)
var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

pca_df_base <- data.frame(
  sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

# --- FUNCIÓN MODIFICADA: Asignar colores basados en repetición de ID ---
# Esta función generará la columna de color para `harmonized_donor_id`
assign_colors_by_repetition <- function(metadata_df, id_column = "harmonized_donor_id") {
  df <- metadata_df %>% arrange(EpiRR) # Aseguramos un orden consistente si es importante
  
  # Obtener IDs únicos y repetidos
  id_counts <- table(df[[id_column]])
  unique_ids <- names(id_counts[id_counts == 1])
  repeated_ids <- names(id_counts[id_counts > 1])
  
  # Colores para IDs repetidos (asegurando suficientes colores)
  num_repeated_ids <- length(repeated_ids)
  # Usar una paleta más amplia para asegurar suficiente contraste
  if (num_repeated_ids > 0) {
    # Genera una paleta de colores diversa para los IDs repetidos
    # Considera usar un colorRampPalette si hay muchos IDs repetidos
    if (num_repeated_ids <= 8) {
      repeated_colors_base <- RColorBrewer::brewer.pal(n = max(3, num_repeated_ids), name = "Set1")
    } else {
      repeated_colors_base <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(num_repeated_ids)
    }
  } else {
    repeated_colors_base <- character(0) # Ningún ID repetido
  }
  
  # Mapeo de IDs repetidos a colores específicos
  color_map <- setNames(repeated_colors_base, repeated_ids)
  
  # Asignar colores
  colors_assigned <- sapply(df[[id_column]], function(id) {
    if (id %in% unique_ids) {
      return("lightgray") # Gris para IDs únicos
    } else {
      return(color_map[id]) # Color específico para IDs repetidos
    }
  }, USE.NAMES = FALSE)
  
  # Crear la nueva columna con los colores
  df$assigned_color <- colors_assigned
  
  # Preparar los niveles del factor para la leyenda si es necesario
  # Queremos que 'lightgray' sea el primero, seguido de los colores asignados
  color_levels <- c("lightgray", repeated_colors_base)
  
  # Si la columna de ID es un factor, queremos que los niveles de la leyenda
  # se muestren de forma significativa.
  # La etiqueta de la leyenda no será el ID directamente, sino el ID junto con su color.
  # No crearemos un factor aquí, dejaremos que ggplot lo maneje como una columna de caracteres,
  # y manejaremos el orden de los colores en scale_color_manual.
  
  return(df)
}

# --- APLICAR LA NUEVA FUNCIÓN A LOS METADATOS FILTRADOS ---
message("Asignando colores a los IDs de donante basados en su repetición...")
monocyte_meta_filtered <- assign_colors_by_repetition(monocyte_meta_filtered, id_column = "harmonized_donor_id")

# --- 5. Función para generar un plot de PCA (SIN LEYENDA en el plot) ---
# Ahora la leyenda se generará aparte para componerla.
generate_pca_plot_core <- function(pca_data_frame, metadata_df, color_column,
                                   title_prefix = "",
                                   point_size = 3, label_points = FALSE,
                                   colors_palette = NULL, # Ahora se usará para los colores calculados
                                   show_legend_for_extraction = FALSE) {
  
  pca_df_final <- pca_data_frame %>%
    left_join(metadata_df %>% select(EpiRR, !!sym(color_column), harmonized_donor_id), # Incluir harmonized_donor_id para la leyenda
              by = c("sample" = "EpiRR")) %>%
    rename(color_group = !!sym(color_column))
  
  pc1_var_exp <- round(100*var_exp[1], 1)
  pc2_var_exp <- round(100*var_exp[2], 1)
  
  # Crear un factor para 'color_group' para que ggplot sepa el orden
  # Los niveles serán los colores mismos
  # Queremos la leyenda mostrando los harmonized_donor_id, no los colores directamente
  
  # Generar las etiquetas para la leyenda: "ID del Donante [Color]"
  # Crear un dataframe temporal para manejar las etiquetas de la leyenda
  legend_data <- pca_df_final %>%
    distinct(harmonized_donor_id, color_group) %>%
    mutate(legend_label = paste0(harmonized_donor_id)) %>% # Las etiquetas serán los IDs de donante
    arrange(desc(color_group == "lightgray")) # Pone gris al principio
  
  # Asegurar que 'color_group' en pca_df_final tenga los niveles correctos para el orden de la leyenda
  pca_df_final$color_group <- factor(pca_df_final$color_group, levels = unique(legend_data$color_group))
  
  # Mapear los colores a las etiquetas para scale_color_manual
  color_values_map <- setNames(legend_data$color_group, legend_data$legend_label)
  
  p <- ggplot(pca_df_final, aes(x = PC1, y = PC2, color = harmonized_donor_id)) + # Colorear por harmonized_donor_id
    geom_point(size = point_size) +
    labs(
      title = paste0(title_prefix, "PCA por Donante (Repetición)"),
      x = paste0("PC1 (", pc1_var_exp, "%)"),
      y = paste0("PC2 (", pc2_var_exp, "%)"),
      color = "Donante ID" # Título de la leyenda
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
  
  # Aplicar los colores manualmente
  # Usamos pca_df_final$color_group para obtener el color real de cada punto
  # y pca_df_final$harmonized_donor_id para la etiqueta en la leyenda.
  # Esto es un poco truculento porque ggplot asume que el factor de 'color' es lo que quieres en la leyenda.
  # Necesitamos un vector con los colores y los nombres de los donantes.
  
  # Aquí pasamos los colores directamente. Si quieres que la leyenda muestre los donantes,
  # y los puntos se colorean según 'assigned_color', necesitas un mapping explícito.
  # La mejor forma es que 'color_group' en el aes sea el color ya asignado
  # y luego controlar las etiquetas de la leyenda con 'labels' en scale_color_manual.
  
  # Para asegurar que la leyenda muestra los IDs de donante y los colores correspondientes:
  # Primero creamos un mapeo de donante ID a su color
  donor_color_mapping <- pca_df_final %>%
    distinct(harmonized_donor_id, color_group) %>%
    # Aseguramos el orden: gris primero
    arrange(desc(color_group == "lightgray"), harmonized_donor_id)
  
  # Creamos el vector de colores con nombres de donantes
  # Filtrar solo donantes que NO son gris
  custom_colors <- donor_color_mapping %>%
    filter(color_group != "lightgray") %>%
    { setNames(.$color_group, .$harmonized_donor_id) }
  
  # Aplicar la escala de color: 
  #  - muestra solo repetidos en la leyenda
  #  - los grises siguen en el gráfico pero no aparecen en la leyenda
  p <- p + scale_color_manual(values = custom_colors, na.value = "lightgray")
  
  return(p)
}


# --- 6. Definir las columnas del metadata por las que quieres iterar ---
# Ahora solo necesitamos procesar "harmonized_donor_id" una vez,
# ya que la lógica de colores ya está precalculada en 'assigned_color'
# y la pasaremos como la columna de color.
message("Identificando columnas adecuadas para colorear los plots...")

# columns_to_plot_by ahora solo contendrá "assigned_color"
columns_to_plot_by <- "assigned_color" # Esta será la columna que contenga los colores directos

# Esta comprobación ya no es tan crítica porque 'assigned_color' siempre existirá
if (length(columns_to_plot_by) == 0) {
  warning("No se encontraron columnas categóricas adecuadas en los metadatos de monocitos para generar plots.
          Asegúrate de tener columnas con al menos dos categorías y no demasiadas.")
} else {
  message(paste("Se generarán plots para la columna:", paste(columns_to_plot_by, collapse = ", ")))
}

# --- 7. Generar y guardar plots con leyenda separada debajo en un único PDF (tamaño A4) ---
composed_plots_list <- list() # Lista para almacenar los plots compuestos
output_dir <- "PCA_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
output_pdf_path <- file.path(output_dir, "PCA_Macrophage_Donor_Repetition_Plots_A4.pdf")

# Las paletas de colores ya se manejan dentro de assign_colors_by_repetition

if (length(columns_to_plot_by) > 0) {
  a4_width_inches <- 8.27
  a4_height_inches <- 11.69
  
  message(paste("Guardando plots en:", output_pdf_path, "con tamaño A4."))
  
  pdf(output_pdf_path, width = a4_width_inches, height = a4_height_inches)
  
  # La iteración se simplifica porque solo tenemos una "columna" de color ahora
  # que es `assigned_color`.
  col_name <- "assigned_color" # La columna que contiene los colores directos
  message(paste("Procesando columna para plot:", col_name))
  
  # generate_pca_plot_core ya no necesita una `colors_palette` externa
  # porque los colores ya están definidos en la columna `assigned_color`
  # y se manejan en scale_color_manual.
  
  # 1. Generar el plot que usaremos para extraer la leyenda
  plot_for_legend_extraction <- generate_pca_plot_core(
    pca_data_frame = pca_df_base,
    metadata_df = monocyte_meta_filtered,
    color_column = col_name, # Usamos 'assigned_color' para el mapeo estético
    title_prefix = "",
    point_size = 2.5,
    label_points = FALSE,
    show_legend_for_extraction = TRUE
  )
  
  # 2. Extraer la leyenda
  plot_legend <- get_legend(plot_for_legend_extraction)
  
  # 3. Generar el plot final SIN leyenda
  main_plot <- generate_pca_plot_core(
    pca_data_frame = pca_df_base,
    metadata_df = monocyte_meta_filtered,
    color_column = col_name, # Usamos 'assigned_color' para el mapeo estético
    title_prefix = "",
    point_size = 2.5,
    label_points = FALSE,
    show_legend_for_extraction = FALSE
  )
  
  # 4. Componer el plot principal y la leyenda usando plot_grid
  composed_plot <- plot_grid(
    main_plot,
    plot_legend,
    ncol = 1,
    rel_heights = c(0.8, 0.2),
    align = "v"
  )
  
  # Imprimir el plot compuesto en el PDF
  print(composed_plot)
  
  dev.off() # Cierra el dispositivo PDF
  message("Plots con leyenda individual y tamaño controlado guardados exitosamente en PDF.")
  
} else {
  message("No se generaron plots debido a la falta de columnas válidas.")
}

message("Proceso completo.")