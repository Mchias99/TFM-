# ==============================================================================
# FIGURAS -- TFM RESISTOMA DE M. tuberculosis
# ==============================================================================
# Parte A: las figuras y tablas que aparecen en la memoria.
# Parte B: figuras exploratorias adicionales, no incluidas en el documento
#          final pero que se conservan por si son útiles para la defensa.
# ==============================================================================

library(tidyverse)
library(treemapify)   # necesario para el treemap de la Parte B

# Este script da por hecho que el directorio de trabajo de R ya apunta a la
# carpeta del proyecto (donde están los .csv y el resto de scripts). En
# RStudio, la forma más rápida es Session > Set Working Directory > To
# Source File Location; si abres el proyecto desde un .Rproj en la raíz,
# esto ya queda fijado solo.

df_resumido  <- read_csv("matriz_resistoma_resumida_muestras_TFM.csv", show_col_types = FALSE)
df_detallado <- read_csv("matriz_resistoma_detallada_TFM.csv",         show_col_types = FALSE)

# ------------------------------------------------------------------------
# Orden y paleta comunes a todas las figuras (una única definición: si se
# cambia aquí, cambia en todo el script de forma consistente)
# ------------------------------------------------------------------------
orden_paises  <- c("USA", "China", "India", "South Africa", "United Kingdom")
orden_periodo <- c("2000-2007", "2008-2015", "2016-2024")
orden_drtype  <- c("Susceptible", "Other", "HR-TB", "RR-TB", "MDR-TB", "Pre-XDR-TB")

colores_drtype <- c(
  "Susceptible" = "#4C7A3D",  # verde: sensible
  "Other"       = "#8C8C8C",  # gris: fuera de la escalera OMS (ver Discusión)
  "HR-TB"       = "#E6B84D",
  "RR-TB"       = "#D97B4A",
  "MDR-TB"      = "#C0392B",
  "Pre-XDR-TB"  = "#7B2D26"
)

df_resumido <- df_resumido %>%
  mutate(
    isolation_country = factor(isolation_country, levels = orden_paises),
    periodo           = factor(periodo, levels = orden_periodo),
    drtype            = factor(drtype, levels = orden_drtype),
    es_resistente     = as.integer(drtype != "Susceptible")
  )

tema_tfm <- theme_minimal(base_size = 13)

# ==============================================================================
# PARTE A: FIGURAS Y TABLAS DE LA MEMORIA
# ==============================================================================

# ------------------------------------------------------------------------
# Tabla 3: cohorte por país y período
# ------------------------------------------------------------------------
tabla3 <- df_resumido %>%
  count(isolation_country, periodo) %>%
  pivot_wider(names_from = periodo, values_from = n, values_fill = 0) %>%
  mutate(n_total = rowSums(across(where(is.numeric))))

write_csv(tabla3, "Tabla3_cohorte_pais_periodo.csv")

# ------------------------------------------------------------------------
# Figura: distribución global del perfil de resistencia (donut)
# ------------------------------------------------------------------------
df_resumido %>%
  count(drtype) %>%
  mutate(porcentaje = round(n / sum(n) * 100, 1)) %>%
  ggplot(aes(x = 2, y = porcentaje, fill = drtype)) +
  geom_col(color = "white", linewidth = 0.7) +
  coord_polar(theta = "y") +
  xlim(0.8, 2.5) +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  geom_text(aes(label = ifelse(porcentaje >= 3, paste0(porcentaje, "%"), "")),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  theme_void() +
  labs(title = "Distribución porcentual global del perfil de resistencia")

ggsave("Figura_donut_resistencia_global.png", width = 7, height = 6, dpi = 300)

# ------------------------------------------------------------------------
# Figura 3: drtype por país (% apilado), "Other" justo detrás de "Susceptible"
# ------------------------------------------------------------------------
df_resumido %>%
  count(isolation_country, drtype) %>%
  group_by(isolation_country) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = isolation_country, y = pct, fill = drtype)) +
  geom_col() +
  scale_fill_manual(values = colores_drtype, name = "Tipo de resistencia") +
  labs(x = NULL, y = "% de muestras",
       title = "Distribución del tipo de resistencia por país") +
  tema_tfm

ggsave("Figura3_drtype_por_pais.png", width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# Figura 4: top 10 genes de resistencia
# ------------------------------------------------------------------------
genes_top10 <- df_detallado %>%
  filter(!is.na(drugs), drugs != "-") %>%
  distinct(sra_clean, gene_name) %>%
  count(gene_name, sort = TRUE) %>%
  slice_max(n, n = 10)

ggplot(genes_top10, aes(x = reorder(gene_name, n), y = n)) +
  geom_col(fill = "#F16913") +
  geom_text(aes(label = n), hjust = -0.3) +
  coord_flip() +
  expand_limits(y = max(genes_top10$n) * 1.1) +
  labs(x = NULL, y = "Nº de muestras (de 225) con mutación de resistencia",
       title = "Genes más implicados en la resistencia (top 10)") +
  tema_tfm

ggsave("Figura4_top_genes_resistencia.png", width = 7, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# Figura: nº de muestras y perfil de resistencia año a año (frecuencias
# absolutas) -- ilustra la disparidad de tamaño muestral entre años que se
# discute como limitación del análisis temporal
# ------------------------------------------------------------------------
df_por_anyo <- df_resumido %>% filter(!is.na(collection_year))
totales_anyo <- df_por_anyo %>% count(collection_year, name = "n_total")

ggplot(df_por_anyo, aes(x = collection_year, fill = drtype)) +
  geom_bar(width = 0.8, color = "black", linewidth = 0.2) +
  geom_text(data = totales_anyo, aes(x = collection_year, y = n_total, label = n_total),
            inherit.aes = FALSE, vjust = -0.5, size = 3, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = seq(2000, 2024, by = 2)) +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(title = "Muestras y perfil de resistencia por año de recolección",
       subtitle = "La altura de cada barra refleja la disparidad de tamaño muestral entre años",
       x = "Año de recolección", y = "Nº de muestras") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("Figura_muestras_por_anyo.png", width = 9, height = 5.5, dpi = 300)

# ------------------------------------------------------------------------
# Figura: % de resistencia por año individual, con línea de tendencia
# ------------------------------------------------------------------------
df_resumido %>%
  group_by(collection_year) %>%
  summarise(n = n(), pct = mean(es_resistente) * 100, .groups = "drop") %>%
  ggplot(aes(x = collection_year, y = pct, size = n)) +
  geom_point(alpha = 0.7, color = "#F16913") +
  geom_smooth(method = "lm", se = FALSE, color = "#C0392B", linetype = "dashed") +
  labs(x = "Año de recolección", y = "% de muestras resistentes", size = "n",
       title = "% de resistencia por año") +
  tema_tfm

ggsave("Figura_resistencia_por_anyo.png", width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# Figura: tendencia temporal de la resistencia, por país (forest plot)
# Se recalcula aquí la regresión por país (en vez de depender de que
# Estadisticaa.R se haya ejecutado antes en la misma sesión), para que este
# script funcione de forma autónoma.
# ------------------------------------------------------------------------
tendencia_por_pais <- map_dfr(orden_paises, function(pais) {
  sub <- df_resumido %>%
    filter(isolation_country == pais) %>%
    mutate(year_c = collection_year - mean(collection_year))
  m <- glm(es_resistente ~ year_c, data = sub, family = binomial)
  tibble(
    pais    = pais,
    OR      = exp(coef(m)["year_c"]),
    IC_2.5  = exp(confint(m)["year_c", 1]),
    IC_97.5 = exp(confint(m)["year_c", 2])
  )
})

tendencia_por_pais %>%
  mutate(pais = fct_reorder(pais, OR)) %>%
  ggplot(aes(x = OR, y = pais)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = IC_2.5, xmax = IC_97.5), height = 0.15, color = "#F16913") +
  geom_point(size = 3, color = "#F16913") +
  scale_x_log10() +
  labs(x = "OR por año (escala log)", y = NULL,
       title = "Tendencia temporal de la resistencia, por país",
       subtitle = "OR < 1: resistencia decreciente en el tiempo  |  OR > 1: creciente") +
  tema_tfm

ggsave("Figura_tendencia_temporal_por_pais.png", width = 8, height = 5, dpi = 300)

cat("\n>>> FIGURAS DE LA MEMORIA GENERADAS <<<\n")

# ==============================================================================
# PARTE B: FIGURAS EXPLORATORIAS (no incluidas en la memoria final)
# ==============================================================================
# Se conservan por si resultan útiles para la defensa oral o como material
# suplementario, pero no tienen un número de figura asignado en el documento.

# B.1 Matriz país x linaje x drtype (balloon plot)
df_resumido %>%
  count(isolation_country, main_lineage, drtype) %>%
  ggplot(aes(x = main_lineage, y = isolation_country)) +
  geom_point(aes(size = n, fill = drtype), shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
  scale_size_continuous(range = c(3, 14), name = "Nº muestras") +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(title = "Matriz país x linaje x resistencia", x = "Linaje principal", y = "País de origen") +
  tema_tfm +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# B.2 Treemap jerárquico país > linaje > drtype
df_resumido %>%
  drop_na(isolation_country, main_lineage, drtype) %>%
  count(isolation_country, main_lineage, drtype) %>%
  filter(n > 0) %>%
  ggplot(aes(area = n, fill = drtype, subgroup = isolation_country, subgroup2 = main_lineage)) +
  geom_treemap() +
  geom_treemap_subgroup_border(color = "black", linewidth = 1.2) +
  geom_treemap_subgroup2_border(color = "grey40", linewidth = 0.6) +
  geom_treemap_subgroup_text(place = "topleft", fontface = "bold", color = "black", alpha = 0.85) +
  geom_treemap_text(aes(label = paste0(drtype, "\nn=", n)), color = "black", place = "center", size = 9) +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(title = "Distribución jerárquica: país > linaje > resistencia") +
  theme(legend.position = "bottom")

# B.3 Evolución porcentual por país y período (dot plot facetado)
df_resumido %>%
  count(periodo, isolation_country, drtype) %>%
  group_by(periodo, isolation_country) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = drtype, y = pct, color = drtype)) +
  geom_segment(aes(xend = drtype, yend = 0), linewidth = 0.8) +
  geom_point(size = 3.5) +
  facet_grid(isolation_country ~ periodo) +
  scale_color_manual(values = colores_drtype) +
  coord_flip() +
  labs(title = "Evolución porcentual de resistencia por país y período",
       x = "Perfil de resistencia", y = "Proporción (%)") +
  tema_tfm +
  theme(legend.position = "none")

# B.4 Top 15 genes, coloreado por drtype (version mas detallada de la Figura 4)
df_genes_largo <- df_resumido %>%
  separate_rows(genes_afectados, sep = ";\\s*") %>%
  filter(!is.na(genes_afectados), genes_afectados != "")

top15_genes <- df_genes_largo %>% count(genes_afectados, sort = TRUE) %>% slice_max(n, n = 15)

df_genes_largo %>%
  filter(genes_afectados %in% top15_genes$genes_afectados) %>%
  ggplot(aes(x = fct_infreq(genes_afectados), fill = drtype)) +
  geom_bar(color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = colores_drtype, name = "Perfil DR") +
  labs(title = "Top 15 genes más mutados, por perfil de resistencia",
       x = "Gen", y = "Nº de mutaciones") +
  tema_tfm

cat(">>> FIGURAS EXPLORATORIAS (PARTE B) GENERADAS <<<\n")
