library(tidyverse)
library(DescTools)


# 0. CARGA DE DATOS Y DEFINICIÓN DE VARIABLES

# Se carga la matriz resumida una sola vez y se crean las variables derivadas
df_resistoma_resumido <- read_csv("matriz_resistoma_resumida_muestras_TFM.csv", show_col_types = FALSE) %>%
  mutate(
    es_resistente = as.integer(drtype != "Susceptible"),
    periodo_num   = as.numeric(as.factor(periodo))
  )

cat("\n--- DISTRIBUCIÓN GLOBAL DE RESISTENCIA (0 = Susceptible, 1 = Resistente) ---\n")
print(table(df_resistoma_resumido$es_resistente))


# 1. ASOCIACIÓN GEOGRÁFICA Y TEMPORAL DE RESISTENCIA (drtype)

# 1.1 Por País (Test Exacto de Fisher con simulación de Monte Carlo)
tabla_pais <- table(df_resistoma_resumido$isolation_country, df_resistoma_resumido$drtype)
cat("\n--- TABLA PAÍS vs DRTYPE ---\n")
print(tabla_pais)

test_fisher_pais <- fisher.test(tabla_pais, simulate.p.value = TRUE, B = 10000)
cat("\n--- FISHER PAÍS vs DRTYPE ---\n")
print(test_fisher_pais)

# 1.2 Por Período (Fisher)
tabla_periodo <- table(df_resistoma_resumido$periodo, df_resistoma_resumido$drtype)
cat("\n--- TABLA PERÍODO vs DRTYPE ---\n")
print(tabla_periodo)

test_fisher_periodo <- fisher.test(tabla_periodo, simulate.p.value = TRUE, B = 10000)
cat("\n--- FISHER PERÍODO vs DRTYPE (SIMULACIÓN MONTE CARLO) ---\n")
print(test_fisher_periodo)


# 2. TENDENCIA TEMPORAL (COCHRAN-ARMITAGE)


# 2.1 Test de tendencia de Cochran-Armitage (alternative = "one.sided")
tabla_tendencia <- table(df_resistoma_resumido$periodo, df_resistoma_resumido$es_resistente)
test_tendencia <- CochranArmitageTest(tabla_tendencia, alternative = "one.sided")
cat("\n--- TEST DE TENDENCIA DE COCHRAN-ARMITAGE ---\n")
print(test_tendencia)

# 2.2 Regresión logística de tendencia temporal continua ajustada por País
modelo_tendencia_ajustado <- glm(
  es_resistente ~ periodo_num + isolation_country,
  data = df_resistoma_resumido,
  family = binomial
)
cat("\n--- REGRESIÓN LOGÍSTICA DE TENDENCIA AJUSTADA POR PAÍS ---\n")
print(summary(modelo_tendencia_ajustado))


# 3. REGRESIÓN LOGÍSTICA MULTIVARIANTE (CATEGÓRICA PAÍS + PERÍODO)

modelo_multivariante <- glm(
  es_resistente ~ isolation_country + periodo,
  data = df_resistoma_resumido,
  family = binomial
)
cat("\n--- REGRESIÓN LOGÍSTICA MULTIVARIANTE (PAÍS + PERÍODO) ---\n")
print(summary(modelo_multivariante))

cat("\n--- ODDS RATIO E INTERVALOS DE CONFIANZA (MODELO MULTIVARIANTE) ---\n")
print(exp(coef(modelo_multivariante)))
print(suppressMessages(exp(confint(modelo_multivariante))))


# 4. CORRECCIÓN DINÁMICA DE P-VALORES (BENJAMINI-HOCHBERG / FDR)

p_valores_extraidos <- c(
  Fisher_Pais        = test_fisher_pais$p.value,
  Fisher_Periodo     = test_fisher_periodo$p.value,
  Tendencia_Periodo  = test_tendencia$p.value
)

p_valores_modelo <- summary(modelo_multivariante)$coefficients[-1, 4]

p_valores_globales <- c(p_valores_extraidos, p_valores_modelo)
p_ajustados        <- p.adjust(p_valores_globales, method = "BH")

cat("\n--- CORRECCIÓN DE P-VALORES (FDR - BENJAMINI & HOCHBERG) ---\n")
print(data.frame(
  Prueba_Coeficiente = names(p_valores_globales),
  P_Original         = round(p_valores_globales, 5),
  P_Ajustado_BH      = round(p_ajustados, 5),
  row.names          = NULL
))


# 5. ANÁLISIS TEMPORAL CONTINUO POR AÑO DE RECOLECCIÓN


# 5.1 Chi-cuadrado exacto por año individual (Monte Carlo)
tabla_temporal <- table(df_resistoma_resumido$collection_year, df_resistoma_resumido$drtype)
test_chi_temporal <- chisq.test(tabla_temporal, simulate.p.value = TRUE, B = 10000)
cat("\n--- RESULTADO CHI-CUADRADO (MONTE CARLO) TEMPORAL POR AÑO ---\n")
print(test_chi_temporal)

# 5.2 Regresión Logística Simple por Año
modelo_logistico_ano <- glm(
  es_resistente ~ collection_year,
  data = df_resistoma_resumido,
  family = binomial
)
cat("\n--- REGRESIÓN LOGÍSTICA (AÑO CONTINUO) ---\n")
print(summary(modelo_logistico_ano))

cat("\n--- ODDS RATIO E IC 95% (AÑO CONTINUO) ---\n")
resultados_or_ano <- exp(cbind(OR = coef(modelo_logistico_ano), suppressMessages(confint(modelo_logistico_ano))))
print(resultados_or_ano)


# 6. TENDENCIA TEMPORAL INDEPENDIENTE POR PAÍS

paises <- c("USA", "China", "India", "South Africa", "United Kingdom")

resultados_pais <- map_dfr(paises, function(p) {
  sub <- df_resistoma_resumido %>%
    filter(isolation_country == p) %>%
    mutate(year_c = collection_year - mean(collection_year, na.rm = TRUE))
  
  m <- glm(es_resistente ~ year_c, data = sub, family = binomial)
  
  tibble(
    pais    = p,
    OR      = exp(coef(m)["year_c"]),
    IC_2.5  = suppressMessages(exp(confint(m)["year_c", 1])),
    IC_97.5 = suppressMessages(exp(confint(m)["year_c", 2])),
    p       = summary(m)$coefficients["year_c", 4]
  )
})

cat("\n--- TENDENCIAS POR PAÍS INDIVIDUAL ---\n")
print(resultados_pais, n = Inf)

write_csv(resultados_pais, "Tendencia_temporal_por_pais.csv")