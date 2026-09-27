# ==============================================================================
# ANÁLISIS ESTADÍSTICO -- TFM RESISTOMA DE M. tuberculosis
# ==============================================================================
# Todo el análisis parte de matriz_resistoma_resumida_muestras_TFM.csv (Fase_3.R).
# Estructura de este script, en el mismo orden en que se citan los resultados
# en la memoria:
#   1. Asociación entre drtype y país / período (chi-cuadrado y Fisher)
#   2. Tendencia temporal (Cochran-Armitage, por período y por año)
#   3. Regresión logística ajustada por país y período
#   4. Comparaciones país a país con corrección por comparaciones múltiples
#   5. Efecto del año como variable continua (global y por país)
# ==============================================================================

library(tidyverse)
library(DescTools)   # CochranArmitageTest

# Este script da por hecho que el directorio de trabajo de R ya apunta a la
# carpeta del proyecto (donde están los .csv y el resto de scripts). En
# RStudio, la forma más rápida es Session > Set Working Directory > To
# Source File Location; si abres el proyecto desde un .Rproj en la raíz,
# esto ya queda fijado solo.

df <- read_csv("matriz_resistoma_resumida_muestras_TFM.csv", show_col_types = FALSE)

orden_paises  <- c("USA", "China", "India", "South Africa", "United Kingdom")
orden_periodo <- c("2000-2007", "2008-2015", "2016-2024")

df <- df %>%
  mutate(
    isolation_country = factor(isolation_country, levels = orden_paises),
    periodo           = factor(periodo, levels = orden_periodo),
    es_resistente     = as.integer(drtype != "Susceptible")   # 1 = resistente, 0 = susceptible
  )

# ==============================================================================
# 1. ASOCIACIÓN ENTRE EL PERFIL DE RESISTENCIA (drtype) Y PAÍS / PERÍODO
# ==============================================================================

# 1.1 Por país -- chi-cuadrado + Fisher simulado (hay celdas con esperado <5,
#     así que el Fisher por simulación de Monte Carlo es el test de referencia)
tabla_pais <- table(df$isolation_country, df$drtype)
print(tabla_pais)
chisq.test(tabla_pais)
fisher.test(tabla_pais, simulate.p.value = TRUE, B = 10000)

# 1.2 Por período
tabla_periodo <- table(df$periodo, df$drtype)
print(tabla_periodo)
chisq.test(tabla_periodo)

# ==============================================================================
# 2. TENDENCIA TEMPORAL DE LA RESISTENCIA (COCHRAN-ARMITAGE)
# ==============================================================================

# 2.1 Por período (3 tramos)
tabla_tendencia_periodo <- table(df$periodo, df$es_resistente)
CochranArmitageTest(tabla_tendencia_periodo, alternative = "increasing")

# 2.2 Por año real, sin agrupar en tramos -- comprobación de robustez con la
#     máxima resolución temporal disponible
tabla_tendencia_anyo <- table(df$collection_year, df$es_resistente)
CochranArmitageTest(tabla_tendencia_anyo, alternative = "two.sided")

# 2.3 Chi-cuadrado exacto por año individual (Monte Carlo). Con 25 años y
#     hasta 6 categorías de drtype, muchas celdas quedan con 0-2 muestras: se
#     reporta como exploratorio, no como prueba principal de tendencia.
tabla_anyo_drtype <- table(df$collection_year, df$drtype)
chisq.test(tabla_anyo_drtype, simulate.p.value = TRUE, B = 10000)

# ==============================================================================
# 3. REGRESIÓN LOGÍSTICA AJUSTADA (PAÍS + PERÍODO)
# ==============================================================================
# Referencia: EE. UU. (menor proporción de resistencia observada) y período
# 2000-2007 (el más antiguo). Sin este releveling, R usa el primer nivel por
# orden alfabético (China) como referencia y los odds ratios no coinciden con
# los de la memoria.
df <- df %>%
  mutate(
    isolation_country = relevel(isolation_country, ref = "USA"),
    periodo            = relevel(periodo, ref = "2000-2007")
  )

modelo_pais_periodo <- glm(es_resistente ~ isolation_country + periodo,
                            data = df, family = binomial)
summary(modelo_pais_periodo)

tabla_modelo <- tibble(
  variable = names(coef(modelo_pais_periodo)),
  OR       = exp(coef(modelo_pais_periodo)),
  IC_2.5   = exp(confint(modelo_pais_periodo))[, 1],
  IC_97.5  = exp(confint(modelo_pais_periodo))[, 2],
  p        = summary(modelo_pais_periodo)$coefficients[, 4]
)
print(tabla_modelo, n = Inf)
write_csv(tabla_modelo, "Tabla_regresion_logistica.csv")

# ==============================================================================
# 4. COMPARACIONES PAÍS A PAÍS (FISHER EXACTO) + CORRECCIÓN BH
# ==============================================================================
pares <- combn(orden_paises, 2, simplify = FALSE)

comparaciones_pais <- map_dfr(pares, function(par) {
  sub   <- df %>% filter(isolation_country %in% par)
  tabla <- table(droplevels(sub$isolation_country), sub$es_resistente)
  test  <- fisher.test(tabla)
  tibble(comparacion = paste(par, collapse = " vs "), p_raw = test$p.value)
})

comparaciones_pais <- comparaciones_pais %>%
  mutate(p_BH = p.adjust(p_raw, method = "BH")) %>%
  arrange(p_BH)

print(comparaciones_pais, n = Inf)
write_csv(comparaciones_pais, "Comparaciones_pais_BH.csv")

# ==============================================================================
# 5. EFECTO DEL AÑO COMO VARIABLE CONTINUA
# ==============================================================================

# 5.1 Modelo simple, sin ajustar -- es el que se cita en "Resistencias en el
#     tiempo" (OR=0.975; IC95% 0.934-1.016; p=0.238)
modelo_anyo_simple <- glm(es_resistente ~ collection_year, data = df, family = binomial)
summary(modelo_anyo_simple)
cat("\nOR por año adicional (sin ajustar):\n")
print(exp(coef(modelo_anyo_simple)["collection_year"]))
print(exp(confint(modelo_anyo_simple)["collection_year", ]))

# 5.2 El mismo modelo pero ajustando por país -- este es el que se cita en la
#     Discusión/Conclusiones como comprobación de robustez (OR=0.973 por año;
#     IC95% 0.929-1.019; p=0.246). El ajuste por país cambia ligeramente la
#     estimación respecto al modelo simple de 5.1: son dos resultados
#     distintos y ambos aparecen citados en el documento, no uno redundante
#     con el otro.
df <- df %>% mutate(year_c = collection_year - mean(collection_year))

modelo_anyo_ajustado <- glm(es_resistente ~ isolation_country + year_c,
                             data = df, family = binomial)
summary(modelo_anyo_ajustado)
cat("\nOR por año adicional (ajustado por país):\n")
print(exp(coef(modelo_anyo_ajustado)["year_c"]))
cat("OR equivalente por década:\n")
print(exp(coef(modelo_anyo_ajustado)["year_c"] * 10))
print(exp(confint(modelo_anyo_ajustado)["year_c", ]))

# 5.3 Modelo independiente por país -- revela si la ausencia de tendencia
#     global esconde trayectorias opuestas entre países (ver Discusión)
tendencia_por_pais <- map_dfr(orden_paises, function(pais) {
  sub <- df %>%
    filter(isolation_country == pais) %>%
    mutate(year_c = collection_year - mean(collection_year))

  m <- glm(es_resistente ~ year_c, data = sub, family = binomial)

  tibble(
    pais    = pais,
    OR      = exp(coef(m)["year_c"]),
    IC_2.5  = exp(confint(m)["year_c", 1]),
    IC_97.5 = exp(confint(m)["year_c", 2]),
    p       = summary(m)$coefficients["year_c", 4]
  )
})

print(tendencia_por_pais, n = Inf)
write_csv(tendencia_por_pais, "Tendencia_temporal_por_pais.csv")

cat("\n>>> ANÁLISIS ESTADÍSTICO COMPLETADO <<<\n")
