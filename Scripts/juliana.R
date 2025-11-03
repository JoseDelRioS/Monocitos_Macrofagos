library(readr)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(data.table)

leeract_raw<-read_csv("~/Monocitos_Macrofagos/Data/expected_counts_deseq2.csv")
ayudaact<-fread("~/RScripts/metadata_monomacro_filtered.csv.gz")


macrophage<- ayudaact %>%
  filter(harmonized_donor_id %in% c(
    "S00BS4", "S00C1H", "S00BHQ", "S00FTN", "S001S7", "S0022I",
    "S0018A", "S00DVR", "S007SK", "S00E8W", "S00390", "S006VI",
    "S001MJ", "S00H6O", "S00622", "S00J8C"
  ))

m0<- macrophage %>%
  filter(harmonized_cell_type == "macrophage")

m0m<- macrophage %>%
  filter(harmonized_cell_type == "macrophage") %>%
  filter(harmonized_donor_sex == "male")
m0f<- macrophage %>%
  filter(harmonized_cell_type == "macrophage") %>%
  filter(harmonized_donor_sex == "female")
m1m<- macrophage %>%
  filter(harmonized_cell_type == "alternatively activated macrophage") %>%
filter(harmonized_donor_sex == "male")
m1f<- macrophage %>%
  filter(harmonized_cell_type == "alternatively activated macrophage") %>%
  filter(harmonized_donor_sex == "female")

m2m<- macrophage %>%
  filter(harmonized_cell_type == "inflammatory macrophage") %>%
  filter(harmonized_donor_sex == "male")


m2f<- macrophage %>%
  filter(harmonized_cell_type == "inflammatory macrophage") %>%
  filter(harmonized_donor_sex == "male")


colsmon <- c("id_col",intersect(m0$EpiRR, colnames(leeract_raw)))
monact <- leeract_raw[, ..colsmon]


leeract <- monact %>%
  dplyr::rowwise() %>%
  dplyr::mutate(mean_other_cols = mean(c_across(!id_col), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::select(id_col, mean_other_cols)



names(leeract)[2] <- "expresion"
leeract$id_col<-sub("\\..*$", "", leeract$id_col)
leerint<- read_tsv("~/RScripts/PCHiC_peak_matrix_cutoff5.tsv")
nombres<-read_csv("~/RScripts/human_genes_grch37_p13.csv")

intMon<-leerint[leerint$Mac0>5, ]


# Librerías
library(biomaRt)



ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl_ids <- leeract$id_col

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

gene_map <- gene_map %>%
  left_join(leeract %>% dplyr::select(id_col, expresion), 
            by = c("ensembl_gene_id" = "id_col"))

inputnB <- gene_map %>%
  dplyr::select(expresion,ensembl_gene_id,hgnc_symbol) %>%
  mutate(
    hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "",
                         ensembl_gene_id,
                         hgnc_symbol)
  ) %>%
  dplyr::select(hgnc_symbol,expresion) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)   # elimina duplicados de hgnc_symbol

gestion <- c(inputnB$hgnc_symbol)    # Pone la primera columna como nombres de fila
inputnB <- inputnB[ , 2, drop = FALSE]      # Deja solo la segunda columna
colnames(inputnB) <- "profile"
rownames(inputnB)<-gestion


genes <- nombres
frags <- intMon
# Filtrar fragmentos válidos
frags <- frags[!is.na(frags$oeChr) & !is.na(frags$oeStart) & !is.na(frags$oeEnd), ]
gr_genes <- GRanges(
  seqnames = genes$chromosome_name,
  ranges = IRanges(start = genes$start_position, end = genes$end_position),
  gene_id = genes$ensembl_gene_id,
  hgnc_symbol = genes$hgnc_symbol
)


matched_bait <- frags %>%
  dplyr::select(fragment_id = baitID, hgnc_symbol = baitName) %>%
  separate_rows(hgnc_symbol, sep = ";") %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol)) %>%
  mutate(fragment_type = "bait")


# Paso 2: Asociación oeID, excluyendo fragmentos ya usados como bait
gr_oe <- GRanges(
  seqnames = frags$oeChr,
  ranges = IRanges(start = frags$oeStart, end = frags$oeEnd),
  fragment_id = frags$oeID
)

hits_oe <- findOverlaps(gr_oe, gr_genes)

matched_oe <- data.frame(
  fragment_id = mcols(gr_oe)$fragment_id[queryHits(hits_oe)],
  ensembl_gene_id = mcols(gr_genes)$gene_id[subjectHits(hits_oe)],
  hgnc_symbol = mcols(gr_genes)$hgnc_symbol[subjectHits(hits_oe)],
  fragment_type = "oe"
)

# Procesar matched_bait
matched_bait_clean <- matched_bait %>%
  mutate(
    fragment_id = paste0("frag", fragment_id)
  ) %>%
  dplyr::select(fragment_id, hgnc_symbol)

# Procesar matched_oe
matched_oe_clean <- matched_oe %>%
  mutate(
    fragment_id = paste0("frag", fragment_id),
    hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "",
                         ensembl_gene_id,
                         hgnc_symbol)
  ) %>%
  dplyr::select(fragment_id, hgnc_symbol)

# Combinar después de limpiar
gf_nb <- bind_rows(matched_bait_clean, matched_oe_clean)

colnames(gf_nb) <- c("V1", "V2")
gf_nb <- gf_nb %>% filter(!is.na(V1) & !is.na(V2))

nBff<-data.frame(intMon$baitID,intMon$oeID)

nBff[] <- lapply(nBff, function(x) paste0("frag", x))
colnames(nBff) <- c("V1", "V2")


# Supone que ya tienes los data frames leerint y nombres cargados en tu entorno

# 1. Selecciona las columnas de interés de leerint y nombres
leerint_bait <- intMon[, c("baitChr", "baitID", "baitStart", "baitEnd")]
leerint_oe   <- intMon[, c("oeChr", "oeID", "oeStart", "oeEnd")]
nombres_sel  <- nombres[, c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name")]

# 2. Construye los data frames individuales con nombres estandarizados
# Del bait
df_bait <- data.frame(
  ID   = paste0("frag", leerint_bait$baitID),
  chr  = paste0("chr", leerint_bait$baitChr),
  start = leerint_bait$baitStart,
  end   = leerint_bait$baitEnd,
  stringsAsFactors = FALSE
)

# Del oe
df_oe <- data.frame(
  ID   = paste0("frag", leerint_oe$oeID),
  chr  = paste0("chr", leerint_oe$oeChr),
  start = leerint_oe$oeStart,
  end   = leerint_oe$oeEnd,
  stringsAsFactors = FALSE
)

# De los genes
# Si hgnc_symbol está vacío, usa ensembl_gene_id
gene_id <- ifelse(
  is.na(nombres_sel$hgnc_symbol) | nombres_sel$hgnc_symbol == "",
  nombres_sel$ensembl_gene_id,
  nombres_sel$hgnc_symbol
)

df_genes <- data.frame(
  ID   = gene_id,
  chr  = paste0("chr", nombres_sel$chromosome_name),
  start = nombres_sel$start_position,
  end   = nombres_sel$end_position,
  stringsAsFactors = FALSE
)

# 3. Une todo en un solo data frame
df_final <- rbind(df_bait, df_oe, df_genes)

# 4. Opcional: elimina filas duplicadas
ann_net_nb <- unique(df_final)

