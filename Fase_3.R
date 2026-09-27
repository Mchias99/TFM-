
# SCRIPT 2: INTEGRACIÓN DEL RESISTOMA 


library(tidyverse)
library(janitor)


# Función de limpieza de cadenas 
limpiar_texto <- function(x) {
  x %>%
    iconv(from = "", to = "UTF-8", sub = "") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT", sub = "") %>%
    trimws()
}

cat("1. Comprobando ficheros de entrada...\n")

archivos_necesarios <- c(
  "metadatos_225_muestras_TFM.csv",
  "tbprofiler_results/tbprofiler.txt",
  "tbprofiler_results/tbprofiler.variants.csv"
)
faltantes <- archivos_necesarios[!file.exists(archivos_necesarios)]
if (length(faltantes) > 0) {
  stop("ERROR: Faltan archivos de entrada: ", paste(faltantes, collapse = ", "))
}


# PASO 1: METADATOS DE LA COHORTE (Fase 1)

df_metadatos <- read_csv("metadatos_225_muestras_TFM.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(across(where(is.character), limpiar_texto)) %>%
  mutate(sra_clean = str_extract(run_accession, "(SRR|ERR|DRR)\\d+|(GCA|GCF)_\\d+\\.\\d+")) %>%
  filter(!is.na(sra_clean)) %>%
  distinct(sra_clean, .keep_all = TRUE)


# PASO 2: RESUMEN TB-PROFILER (Lineajes + Tipo de Resistencia)

df_resumen_tb <- read_tsv("tbprofiler_results/tbprofiler.txt", show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(across(where(is.character), limpiar_texto)) %>%
  mutate(sra_clean = str_extract(sample, "(SRR|ERR|DRR)\\d+|(GCA|GCF)_\\d+\\.\\d+")) %>%
  filter(!is.na(sra_clean)) %>%
  distinct(sra_clean, .keep_all = TRUE) %>%
  select(sra_clean, main_lineage, sub_lineage, drtype)


# PASO 3: DETALLE DE VARIANTES

df_tbprofiler_var <- read_csv("tbprofiler_results/tbprofiler.variants.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(across(where(is.character), limpiar_texto)) %>%
  mutate(sra_clean = str_extract(sample, "(SRR|ERR|DRR)\\d+|(GCA|GCF)_\\d+\\.\\d+")) %>%
  filter(!is.na(sra_clean))


# PASO 4: CONSTRUCCIÓN DE LA TABLA MAESTRA (LEFT_JOIN PARA MANTENER N=225)

cat("2. Construyendo la tabla maestra por muestra (incluye pansensibles)...\n")

df_maestro <- df_metadatos %>%
  left_join(df_resumen_tb, by = "sra_clean") %>%
  select(
    sra_clean, study_accession, isolation_country, collection_year,
    periodo, instrument_model, main_lineage, sub_lineage, drtype
  )

cat("   Muestras en la cohorte con metadatos (N):", nrow(df_maestro), "\n")


# PASO 5: UNIÓN CON DETALLE DE VARIANTES

cat("3. Añadiendo el detalle de variantes (left_join preserva pansensibles)...\n")

df_resistoma_detallado <- df_maestro %>%
  left_join(
    df_tbprofiler_var %>% select(sra_clean, gene_name, change, freq, type, drugs),
    by = "sra_clean"
  ) %>%
  distinct(sra_clean, gene_name, change, drugs, .keep_all = TRUE)

write_csv(df_resistoma_detallado, "matriz_resistoma_detallada_TFM.csv")


# PASO 6: RESUMEN POR MUESTRA 

cat("4. Generando resumen por muestra ...\n")

df_resistoma_resumido <- df_resistoma_detallado %>%
  group_by(sra_clean) %>%
  summarise(
    study_accession   = first(study_accession),
    isolation_country = first(isolation_country),
    collection_year   = first(collection_year),
    periodo           = first(periodo),
    instrument_model  = first(instrument_model),
    main_lineage      = first(main_lineage),
    sub_lineage       = first(sub_lineage),
    drtype            = first(drtype),
    
    total_mutaciones     = sum(!is.na(gene_name)),
    genes_afectados      = paste(unique(na.omit(gene_name)), collapse = "; "),
    mutaciones_cambios   = paste(unique(na.omit(change)), collapse = "; "),
    farmacos_resistencia = paste(unique(na.omit(drugs[drugs != "-" & drugs != "" & !is.na(drugs)])), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    farmacos_resistencia = if_else(farmacos_resistencia == "", "Sensible / Sin resistencia", farmacos_resistencia),
    genes_afectados      = if_else(genes_afectados == "", NA_character_, genes_afectados)
  )

write_csv(df_resistoma_resumido, "matriz_resistoma_resumida_muestras_TFM.csv")


# PASO 7: CÁLCULO DE FRECUENCIAS Y PORCENTAJES (NUEVA TABLA)

cat("5. Calculando tabla global porcentual (Linajes, DR_type y Fármacos)...\n")


total_muestras <- nrow(df_resistoma_resumido)


df_pct_linaje <- df_resistoma_resumido %>%
  count(main_lineage, name = "N") %>%
  mutate(Categoria = "Linaje Principal",
         Porcentaje = round((N / total_muestras) * 100, 2)) %>%
  rename(Subcategoria = main_lineage)


df_pct_drtype <- df_resistoma_resumido %>%
  count(drtype, name = "N") %>%
  mutate(Categoria = "Perfil de Resistencia",
         Porcentaje = round((N / total_muestras) * 100, 2)) %>%
  rename(Subcategoria = drtype)


df_pct_farmacos <- df_resistoma_resumido %>%
  filter(farmacos_resistencia != "Sensible / Sin resistencia") %>%
  separate_rows(farmacos_resistencia, sep = ";\\s*") %>%
  count(farmacos_resistencia, name = "N") %>%
  mutate(Categoria = "Fármaco Específico",
         Porcentaje = round((N / total_muestras) * 100, 2)) %>%
  rename(Subcategoria = farmacos_resistencia)


n_sensibles <- sum(df_resistoma_resumido$farmacos_resistencia == "Sensible / Sin resistencia")
df_pct_sensibles <- tibble(
  Subcategoria = "Sensible / Sin resistencia",
  N = n_sensibles,
  Categoria = "Fármaco Específico",
  Porcentaje = round((n_sensibles / total_muestras) * 100, 2)
)

# 7.4 Unir todo en una única tabla final 
df_porcentajes_global <- bind_rows(
  df_pct_linaje,
  df_pct_drtype,
  df_pct_sensibles,
  df_pct_farmacos
) %>%
  select(Categoria, Subcategoria, N, Porcentaje) %>%
  arrange(Categoria, desc(N)) # Ordena de mayor a menor frecuencia dentro de cada categoría

write_csv(df_porcentajes_global, "tabla_frecuencias_porcentajes_TFM.csv")


cat(" INTEGRACIÓN DEL RESISTOMA COMPLETADA CON ÉXITO:\n")
cat(" Muestras en la cohorte final (N):", total_muestras, "\n")
cat(" Muestras pansensibles:", n_sensibles, "(", round((n_sensibles/total_muestras)*100, 1), "%)\n")
cat(" Archivos exportados:\n")
cat("  1. matriz_resistoma_detallada_TFM.csv (Variantes por muestra)\n")
cat("  2. matriz_resistoma_resumida_muestras_TFM.csv (1 fila = 1 muestra)\n")
cat("  3. tabla_frecuencias_porcentajes_TFM.csv (Métricas globales)\n")

View(df_porcentajes_global)