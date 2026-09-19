# ==============================================================================
# SCRIPT 3a: FIGURAS -- RESULTADOS DEL TFM
# ==============================================================================

library(tidyverse)

setwd("C:/Users/migue/Desktop/TFM_R")

df_resistoma_resumido  <- read_csv("matriz_resistoma_resumida_muestras_TFM.csv", show_col_types = FALSE)
df_resistoma_detallado <- read_csv("matriz_resistoma_detallada_TFM.csv",         show_col_types = FALSE)

orden_paises  <- c("USA", "China", "India", "South Africa", "United Kingdom")
orden_drtype  <- c("Susceptible", "Other", "HR-TB", "RR-TB", "MDR-TB", "Pre-XDR-TB")

df_resistoma_resumido <- df_resistoma_resumido %>%
  mutate(
    isolation_country = factor(isolation_country, levels = orden_paises),
    drtype            = factor(drtype, levels = orden_drtype),
    es_resistente     = as.integer(drtype != "Susceptible")
  )

# ------------------------------------------------------------------------
# FIGURA 3: drtype por pais (% apilado) -- "Other" justo detras de "Susceptible"
# ------------------------------------------------------------------------
colores_drtype <- c("Susceptible" = "#4C7A3D", "Other" = "#8C8C8C", "HR-TB" = "#E6B84D",
                    "RR-TB" = "#D97B4A", "MDR-TB" = "#C0392B", "Pre-XDR-TB" = "#7B2D26")

df_resistoma_resumido %>%
  count(isolation_country, drtype) %>%
  group_by(isolation_country) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = isolation_country, y = pct, fill = drtype)) +
  geom_col() +
  scale_fill_manual(values = colores_drtype) +
  labs(x = NULL, y = "% de muestras", fill = "Tipo de resistencia",
       title = "Distribución del tipo de resistencia por país") +
  theme_minimal(base_size = 13)

ggsave("Figura3_drtype_por_pais.png", width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# FIGURA 4: top 10 genes de resistencia
# ------------------------------------------------------------------------
genes_top <- df_resistoma_detallado %>%
  filter(!is.na(drugs), drugs != "-") %>%
  distinct(sra_clean, gene_name) %>%
  count(gene_name, sort = TRUE) %>%
  slice_max(n, n = 10)

ggplot(genes_top, aes(x = reorder(gene_name, n), y = n)) +
  geom_col(fill = "#F16913") +
  geom_text(aes(label = n), hjust = -0.3) +
  coord_flip() +
  expand_limits(y = max(genes_top$n) * 1.1) +
  labs(x = NULL, y = "Nº de muestras (de 225) con mutación de resistencia",
       title = "Genes más implicados en la resistencia (top 10)") +
  theme_minimal(base_size = 13)

ggsave("Figura4_top_genes_resistencia.png", width = 7, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# FIGURA (exploratoria): % de resistencia por año individual
# ------------------------------------------------------------------------
df_resistoma_resumido %>%
  group_by(collection_year) %>%
  summarise(n = n(), pct = mean(es_resistente) * 100, .groups = "drop") %>%
  ggplot(aes(x = collection_year, y = pct, size = n)) +
  geom_point(alpha = 0.7, color = "#F16913") +
  geom_smooth(method = "lm", se = FALSE, color = "#C0392B", linetype = "dashed", show.legend = FALSE) +
  labs(x = "Año de recolección", y = "% de muestras resistentes", size = "n",
       title = "% de resistencia por año") +
  theme_minimal(base_size = 13)

ggsave("Fig_resistencia_por_anyo.png", width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# Figura: forest plot de la tendencia anual por pais
# ------------------------------------------------------------------------
resultados_pais %>%
  mutate(pais = fct_reorder(pais, OR)) %>%
  ggplot(aes(x = OR, y = pais)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = IC_2.5, xmax = IC_97.5), height = 0.15, color = "#F16913") +
  geom_point(size = 3, color = "#F16913") +
  scale_x_log10() +
  labs(x = "OR por año (escala log)", y = NULL,
       title = "Tendencia temporal de la resistencia, por país",
       subtitle = "OR < 1: resistencia decreciente en el tiempo  |  OR > 1: creciente") +
  theme_minimal(base_size = 13)

ggsave("Figura_tendencia_temporal_por_pais.png", width = 8, height = 5, dpi = 300)


# 1. Nueva paleta de alto contraste (fácil de distinguir)
colores_claros <- c(
  "Susceptible" = "#4DAF4A",  # Verde (Sensible)
  "Other"       = "#377EB8",  # Azul
  "HR-TB"       = "#FF7F00",  # Naranja
  "RR-TB"       = "#FFFF33",  # Amarillo
  "MDR-TB"      = "#E41A1C",  # Rojo (Peligro)
  "Pre-XDR-TB"  = "#984EA3",  # Morado (Alerta máxima)
  "XDR-TB"      = "#A65628"   # Marrón
)

# 2. Filtrado de datos brutos (Rango 2000-2022)
df_figura_absoluta <- df_resistoma_resumido %>%
  filter(!is.na(collection_year) & collection_year >= 2000 & collection_year <= 2022)

# 3. Calcular los totales por año para poner la etiqueta arriba
df_totales_anuales <- df_figura_absoluta %>%
  count(collection_year, name = "n_total")

# 4. Generación del gráfico de frecuencias absolutas
figura_disparidad <- ggplot(df_figura_absoluta, aes(x = collection_year, fill = drtype)) +
  geom_bar(width = 0.8, color = "black", linewidth = 0.2) +
  
  # Etiqueta con el 'n' total sobre cada barra
  geom_text(
    data = df_totales_anuales,
    aes(x = collection_year, y = n_total, label = paste0("n=", n_total)),
    inherit.aes = FALSE,
    vjust = -0.5, size = 3, fontface = "bold"
  ) +
  
  # Escalas y colores
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Da espacio arriba para el texto
  scale_x_continuous(breaks = seq(2000, 2022, by = 2)) +
  scale_fill_manual(values = colores_claros, name = "Perfil DR") +
  
  # Títulos
  labs(
    title = "Distribución Temporal del Resistoma (Frecuencias Absolutas)",
    subtitle = "La altura de las barras refleja la marcada disparidad en el tamaño muestral anual",
    x = "Año de recolección (collection_year)",
    y = "Número total de aislados (n)"
  ) +
  
  # Tema visual limpio
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(figura_disparidad)

# ==============================================================================
# BLOQUE 2: ESTRUCTURACIÓN DE FACTORES, PALETA Y TEMA VISUAL
# ==============================================================================

orden_paises   <- c("USA", "China", "India", "South Africa", "United Kingdom")
orden_periodo  <- c("2000-2007", "2008-2015", "2016-2024")
orden_drtype   <- c("Susceptible", "Other", "HR-TB", "RR-TB", "MDR-TB", "Pre-XDR-TB", "XDR-TB")

colores_drtype <- c(
  "Susceptible" = "#FFF5EB",
  "Other"       = "#FDD0A2",
  "HR-TB"       = "#FDAE6B",
  "RR-TB"       = "#F16913",
  "MDR-TB"      = "#D94801",
  "Pre-XDR-TB"  = "#A63603",
  "XDR-TB"      = "#7F2704"
)

df_resistoma_resumido <- df_resistoma_resumido %>%
  mutate(
    isolation_country = factor(isolation_country, levels = orden_paises),
    periodo           = factor(periodo, levels = orden_periodo),
    drtype            = factor(drtype, levels = intersect(orden_drtype, unique(drtype))),
    es_resistente     = as.integer(drtype != "Susceptible")
  )

tema_tfm <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "#7F2704"),
    axis.title = element_text(face = "bold", size = 10, color = "#7F2704"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold", size = 10, color = "#7F2704"),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# BLOQUE 3: TABLA ESTADÍSTICA DE LA COHORTE
# ==============================================================================

# TABLA 3: Distribución por País y Período
tabla3 <- df_resistoma_resumido %>%
  count(isolation_country, periodo) %>%
  pivot_wider(names_from = periodo, values_from = n, values_fill = 0) %>%
  mutate(n_total = rowSums(across(where(is.numeric))))

write_csv(tabla3, "Tabla3_cohorte_pais_periodo.csv")

# ==============================================================================
# BLOQUE 4: VISUALIZACIÓN DE RESULTADOS (FIGURAS PARA MEMORIA)
# ==============================================================================

# 4.1 FIGURA MACRO: DONUT CHART (GRÁFICO DE RUEDA DE RESISTENCIA GLOBAL)
df_rueda <- df_resistoma_resumido %>%
  count(drtype) %>%
  mutate(porcentaje = round((n / sum(n)) * 100, 1))

ggplot(df_rueda, aes(x = 2, y = porcentaje, fill = drtype)) +
  geom_col(color = "white", linewidth = 0.7) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.8, 2.5) +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  geom_text(
    aes(label = ifelse(porcentaje >= 3, paste0(porcentaje, "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold"
  ) +
  theme_void() +
  labs(title = "Distribución Porcentual Global del Perfil de Resistencia") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "#7F2704"),
    legend.title = element_text(face = "bold", color = "#7F2704")
  )

# 4.2 FIGURA MESO 1: MATRIZ ESPACIO-TEMPORAL (BALLOON PLOT)
df_resistoma_resumido %>%
  count(isolation_country, main_lineage, drtype) %>%
  ggplot(aes(x = main_lineage, y = isolation_country)) +
  geom_point(aes(size = n, fill = drtype), shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
  scale_size_continuous(range = c(3, 14), name = "Nº Muestras") +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(
    title = "Matriz Espacio-Temporal: País vs. Linaje y Resistencia",
    x = "Linaje Principal", y = "País de Origen"
  ) +
  tema_tfm +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  )

# 4.3 FIGURA MESO 2: MAPA JERÁRQUICO (TREEMAP CORREGIDO SIN AREA=0)
df_treemap <- df_resistoma_resumido %>%
  drop_na(isolation_country, main_lineage, drtype) %>%
  count(isolation_country, main_lineage, drtype) %>%
  filter(n > 0)

ggplot(df_treemap, aes(area = n, fill = drtype, subgroup = isolation_country, subgroup2 = main_lineage)) +
  geom_treemap() +
  geom_treemap_subgroup_border(color = "black", size = 1.2) +
  geom_treemap_subgroup2_border(color = "grey40", size = 0.6) +
  geom_treemap_subgroup_text(place = "topleft", fontface = "bold", color = "black", alpha = 0.85, grow = FALSE) +
  geom_treemap_text(aes(label = paste(drtype, paste0("n=", n), sep = "\n")), color = "black", place = "center", size = 9) +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(title = "Distribución Proporcional Jerárquica: País > Linaje > Resistencia") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "#7F2704"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = "#7F2704")
  )

# 4.4 FIGURA MESO 3: EVOLUCIÓN PORCENTUAL ESPACIO-TEMPORAL (DOT PLOT FACETADO)
df_resistoma_resumido %>%
  count(periodo, isolation_country, drtype) %>%
  group_by(periodo, isolation_country) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = drtype, y = pct, color = drtype)) +
  geom_segment(aes(xend = drtype, yend = 0), size = 0.8) +
  geom_point(size = 3.5) +
  facet_grid(isolation_country ~ periodo) +
  scale_color_manual(values = colores_drtype) +
  coord_flip() +
  labs(
    title = "Evolución Porcentual de Resistencia por País y Período",
    x = "Perfil de Resistencia", y = "Proporción (%)"
  ) +
  tema_tfm +
  theme(legend.position = "none")

# 4.5 FIGURA MICRO: TOP 15 GENES MÁS MUTADOS EN EL RESISTOMA
df_genes <- df_resistoma_resumido %>%
  separate_rows(genes_afectados, sep = ";\\s*") %>%
  filter(!is.na(genes_afectados) & genes_afectados != "")

top_genes <- df_genes %>%
  count(genes_afectados, sort = TRUE) %>%
  slice_max(n, n = 15)

ggplot(filter(df_genes, genes_afectados %in% top_genes$genes_afectados), 
       aes(x = fct_infreq(genes_afectados), fill = drtype)) +
  geom_bar(color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(
    title = "Top 15 Genes Más Frecuentemente Mutados en el Resistoma",
    x = "Gen Afectado", y = "Conteo de Mutaciones"
  ) +
  tema_tfm

cat("\n>>> PROCESO COMPLETADO EXITOSAMENTE <<<\n")