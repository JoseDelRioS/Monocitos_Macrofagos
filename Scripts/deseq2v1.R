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
library(DESeq2)

# --- 1. Carga de datos ---
message("Cargando datos...")
leeract_raw <- fread("~/RScripts/genes_tpm_monomacro_processed.csv.gz")
ayudaact <- fread("~/RScripts/metadata_monomacro_filtered.csv.gz")
DESeq2_data <- fread("~/RScripts/genes_expected_count_DESeq2.csv.gz")

monocyte_meta <- ayudaact %>%
  filter(harmonized_donor_id %in% c(
    "S00BS4", "S00C1H", "S00BHQ", "S00FTN", "S001S7", "S0022I",
    "S0018A", "S00DVR", "S007SK", "S00E8W", "S00390", "S006VI",
    "S001MJ", "S00H6O", "S00622", "S00J8C"
  )) %>%
  mutate(across(where(is.character), as.factor))

# --- 3. Preparación de la matriz de expresión ---
message("Preparando matriz de expresión...")
common_samples <- intersect(monocyte_meta$EpiRR, colnames(leeract_raw))

expr_matrix_raw <- leeract_raw %>%
  dplyr::select(id_col, all_of(common_samples))

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


# --- 6. Preparamos los datos ---
datos_meta_DESeq2 <- monocyte_meta_filtered %>%
  dplyr::select(EpiRR, harmonized_cell_type, harmonized_tissue_type, harmonized_donor_sex) %>%
  mutate(
    across(where(is.character), as.factor),
    across(where(is.factor), ~ fct_explicit_na(., na_level = "Missing")),
  )

# --- 6. Preparamos los datos (Continuación) ---

# (Tu código hasta la creación de conteos_DESeq2 es correcto)

# Tu código original:
limpiar_nombres <- function(nombres) {
  nombres_acortados <- gsub("^([^.]*\\.[^.]*)\\..*", "\\1", nombres)
  nombres_unicos <- make.unique(nombres_acortados, sep = "_")
  return(nombres_unicos)
}

conteos_DESeq2_raw <- DESeq2_data %>%
  rename_with(limpiar_nombres) 

# Paso 1: Encontrar las muestras comunes entre los metadatos y los conteos.
# ¡Esto ya lo tenías bien!
common_samples_final <- intersect(datos_meta_DESeq2$EpiRR, colnames(conteos_DESeq2_raw))

# Paso 2: Filtrar AMBOS data frames para que contengan SÓLO estas muestras comunes.
# Y asegurarnos de que estén en el mismo orden.

# Filtra los metadatos
meta_final <- datos_meta_DESeq2 %>%
  filter(EpiRR %in% common_samples_final) %>%
  arrange(match(EpiRR, common_samples_final)) # Ordenar

# Filtra los conteos
conteos_final <- conteos_DESeq2_raw %>%
  dplyr::select(id_col, all_of(common_samples_final)) # Seleccionar y ordenar implícitamente

# --- 7. Convertir a los formatos correctos que DESeq2 exige ---

# Para los METADATOS (colData):
# Los rownames deben ser los identificadores de las muestras.
colData_final <- as.data.frame(meta_final)
rownames(colData_final) <- colData_final$EpiRR

# Para los CONTEOS (countData):
# Debe ser una matriz numérica, con los nombres de genes como rownames.
countData_final <- as.data.frame(conteos_final)
rownames(countData_final) <- countData_final$id_col
countData_final$id_col <- NULL # Eliminar la columna de genes que ahora son rownames
# ¡Importante! Asegúrate de que las cuentas son enteros.
# DESeq2_data a veces se lee con números decimales (ej. 150.0). round() lo soluciona.
countData_final <- round(as.matrix(countData_final)) 

# --- 8. Comprobación Final de Sanidad (El paso que detecta este error) ---
message("Verificando la consistencia final de los datos...")

print(paste("Columnas en countData:", ncol(countData_final)))
print(paste("Filas en colData:", nrow(colData_final)))

if (!all(colnames(countData_final) == rownames(colData_final))) {
  stop("¡ERROR FATAL! Los nombres de las columnas de los conteos no coinciden con los de las filas de los metadatos.")
} else {
  message("¡Éxito! Los datos están alineados y listos para DESeq2.")
}

# --- 9. Crear el objeto DESeqDataSet (Forma Correcta) ---

dds <- DESeqDataSetFromMatrix(countData = countData_final,
                              colData = colData_final,
                              design = ~ harmonized_donor_sex + harmonized_cell_type + harmonized_tissue_type) # El diseño usa la variable 'PC' de los metadatos

message("Objeto DESeqDataSet creado exitosamente.")


# --- 10. Ejecutar el análisis ---
message("Ejecutando el análisis DESeq2...")
dds <- DESeq(dds)
message("Análisis DESeq2 completado.")

# --- 11. Resultados ---
m2vsmo<-results(dds, contrast = c("harmonized_cell_type","alternatively activated macrophage","macrophage"), alpha = 0.05)
m1vsmo<-results(dds, contrast = c("harmonized_cell_type","inflammatory macrophage","macrophage"), alpha = 0.05)



message("Filtrando resultados para el diagrama de Venn...")

# Instalar ggVennDiagram si no está instalado
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}
library(ggVennDiagram)
library(ggplot2)

# Filtrar los resultados de m2vsmo para padj > 0.05 y eliminar filas con NA
genes_m2vsmo_no_significativos <- subset(m2vsmo, padj < 0.05 & !is.na(padj))
# Obtener los nombres de los genes (que están como rownames)
genes_lista_m2 <- rownames(genes_m2vsmo_no_significativos)

# Filtrar los resultados de m1vsmo para padj > 0.05 y eliminar filas con NA
genes_m1vsmo_no_significativos <- subset(m1vsmo, padj < 0.05 & !is.na(padj))
# Obtener los nombres de los genes
genes_lista_m1 <- rownames(genes_m1vsmo_no_significativos)

# --- 13. Creación del diagrama de Venn ---
message("Creando el diagrama de Venn...")

# Crear una lista con los dos conjuntos de genes
lista_genes_venn <- list(
  "M2vsM0 (padj > 0.05)" = genes_lista_m2,
  "M1vsM0 (padj > 0.05)" = genes_lista_m1
)

# Generar el diagrama de Venn
venn_plot <- ggVennDiagram(
  lista_genes_venn,
  set_size = 2.5,
  label_alpha = 0, # Elimina el fondo de los números
  category.names = c("M2vsM0", "M1vsM0") # Nombres más cortos para las categorías
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + # Esquema de color
  labs(title = "Comparación entre tipos de macrófagos vs macrófagos base",
       subtitle = "padj > 0.05")

# Mostrar el gráfico
print(venn_plot)

message("Diagrama de Venn generado exitosamente.")


vsd <- vst(dds, blind = FALSE)
variables_a_graficar <- c("harmonized_cell_type", "harmonized_tissue_type", "harmonized_donor_sex")

pdf(file.path("PCA_Plots", "PCAnuevos.pdf"), width = 11, height = 8.5) # Dimensiones para una página apaisada estándar


for (variable_actual in variables_a_graficar) {

  message(paste("Generando PCA para la variable:", variable_actual))
  pcaData <- plotPCA(vsd, intgroup = variable_actual, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[variable_actual]])) +
    geom_point(size = 4) +
    labs(
      title = "Análisis de Componentes Principales (PCA)",
      subtitle = paste("Muestras coloreadas por:", gsub("_", " ", variable_actual)), # Subtítulo dinámico y legible
      x = paste0("PC1: ", percentVar[1], "% de varianza"),
      y = paste0("PC2: ", percentVar[2], "% de varianza"),
      color = gsub("_", " ", variable_actual) # Título de la leyenda dinámico
    ) +
    theme_bw(base_size = 14) + # Aumentar el tamaño base de la fuente
    coord_fixed()

  # Imprimir el gráfico. Esto es crucial para que se dibuje en el archivo PDF.
  print(pca_plot)

} # Fin del bucle

# 4. Cerrar el dispositivo PDF. Esto guarda el archivo en tu directorio de trabajo.
dev.off()

message("Archivo 'PCAsinPC.pdf' creado exitosamente con todos los gráficos.")

