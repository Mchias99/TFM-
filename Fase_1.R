
# Fase 1: DESCARGA, FILTRADO Y MUESTREO DE LA COHORTE (N = 225)


library(tidyverse)

set.seed(123)

# Este script escribe sus ficheros de salida (accesiones_sra.txt,
# metadatos_225_muestras_TFM.csv) en el directorio de trabajo actual de R.
# Fija antes ese directorio a la carpeta del proyecto -- en RStudio, Session >
# Set Working Directory > To Source File Location -- para que Fase_2.sh y
# Fase_3.R encuentren esos ficheros después sin tener que moverlos a mano.


# PASO 1: DESCARGA DE METADATOS DESDE ENA

cat("1. Descargando metadatos completos con campo de Estudio (study_accession)...\n")

url_ena <- "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_tree(1773)&fields=run_accession,study_accession,country,collection_date,instrument_model,library_strategy,read_count,fastq_ftp&limit=0&format=tsv"

metadatos_raw <- read_tsv(url_ena, show_col_types = FALSE)


# PASO 2: LIMPIEZA, NORMALIZACION Y FILTRADO 

cat("2. Procesando y filtrando metadatos...\n")

metadatos_limpios <- metadatos_raw %>%
  mutate(
    collection_year = as.numeric(str_extract(collection_date, "\\b(19|20)\\d{2}\\b")),
    country_clean = str_trim(str_remove(country, ":.*"))
  ) %>%
  filter(
    collection_year >= 2000 & collection_year <= 2024,
    str_detect(instrument_model, "(?i)Illumina|HiSeq|MiSeq|NextSeq|NovaSeq"),
    str_detect(library_strategy, "(?i)WGS|genomic"),
    as.numeric(read_count) >= 500000,
    !is.na(fastq_ftp) & fastq_ftp != "",
    # Excluir explicitamente valores no válidos o faltantes en país
    !is.na(country_clean),
    country_clean != "",
    !str_detect(country_clean, "(?i)missing|not collected|not provided|unknown")
  ) %>%
  # Normalización de sinonimias de países
  mutate(
    isolation_country = case_when(
      str_detect(country_clean, "(?i)^USA$|United States") ~ "USA",
      str_detect(country_clean, "(?i)^UK$|United Kingdom|Great Britain") ~ "United Kingdom",
      str_detect(country_clean, "(?i)Russia|Russian Federation") ~ "Russia",
      str_detect(country_clean, "(?i)^South Korea$|Republic of Korea") ~ "South Korea",
      str_detect(country_clean, "(?i)^North Korea$|Democratic People.?s Republic of Korea") ~ "North Korea",
      str_detect(country_clean, "(?i)Viet") ~ "Vietnam",
      str_detect(country_clean, "(?i)Democratic Republic of.*Congo") ~ "DR Congo",
      str_detect(country_clean, "(?i)^Congo$|Republic of.*Congo") & !str_detect(country_clean, "(?i)Democratic") ~ "Republic of Congo",
      TRUE ~ country_clean
    )
  ) %>%
  mutate(
    periodo = case_when(
      collection_year >= 2000 & collection_year <= 2007 ~ "2000-2007",
      collection_year >= 2008 & collection_year <= 2015 ~ "2008-2015",
      collection_year >= 2016 & collection_year <= 2024 ~ "2016-2024"
    )
  )

cat("\n--- CONTEO TOP 10 PAÍSES TRAS LIMPIEZA DE SUBREGIONES ---\n")
print(head(sort(table(metadatos_limpios$isolation_country), decreasing = TRUE),10))


# PASO 3: DEDUPLICACIÓN POR PROYECTO

cat("\n3. Controlando sesgo de proyecto (máximo 3 muestras por estudio)...\n")

metadatos_heterogeneos <- metadatos_limpios %>%
  group_by(isolation_country, periodo, study_accession) %>%
  slice_sample(n = 3, replace = FALSE) %>%   
  ungroup()


# PASO 4: SELECCIÓN DE PAÍSES REPRESENTATIVOS (TOP 5 CON MAYOR COBERTURA)


paises_conteo <- metadatos_heterogeneos %>%
  group_by(isolation_country, periodo) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = periodo, values_from = n, values_fill = 0) %>%
  filter(`2000-2007` >= 15 & `2008-2015` >= 15 & `2016-2024` >= 15)

top_paises <- paises_conteo %>%
  mutate(total = `2000-2007` + `2008-2015` + `2016-2024`) %>%
  arrange(desc(total)) %>%
  slice_head(n = 5) %>%
  pull(isolation_country)

cat("\n--- PAÍSES SELECCIONADOS PARA ALCANZAR N = 225 ---\n")
print(top_paises)


# PASO 5: MUESTREO FINAL EXACTO (5 países x 3 períodos x 15 muestras = 225)

dataset_225_definitivo <- metadatos_heterogeneos %>%
  filter(isolation_country %in% top_paises) %>%
  group_by(isolation_country, periodo) %>%
  slice_sample(n = 15, replace = FALSE) %>%
  ungroup()

cat("\n--- COHORTE FINAL DEFINITIVA ---\n")
cat("Total de muestras seleccionadas:", nrow(dataset_225_definitivo), "\n")


# EXPORTACIÓN DE ARCHIVOS

# 1. Lista de accesos SRA para el pipeline de descarga 
write_lines(dataset_225_definitivo$run_accession, "accesiones_sra.txt")

# 2. Metadatos completos
write_csv(dataset_225_definitivo, "metadatos_225_muestras_TFM.csv")

cat("\nArchivos exportados exitosamente con N =", nrow(dataset_225_definitivo), "\n")
cat(" - accesiones_sra.txt (para Fase_2.sh)\n")
cat(" - metadatos_225_muestras_TFM.csv (para Fase_3.R)\n")
