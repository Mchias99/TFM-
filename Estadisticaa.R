# 1. Chi-cuadrado / Fisher -- comparar drtype entre paises o periodos
tabla <- table(df_resistoma_resumido$isolation_country, df_resistoma_resumido$drtype)
print(tabla)
chisq.test(tabla)
fisher.test(tabla, simulate.p.value = TRUE, B = 10000)

# Lo mismo pero por periodo en vez de pais:
tabla_periodo <- table(df_resistoma_resumido$periodo, df_resistoma_resumido$drtype)
print(tabla_periodo)
chisq.test(tabla_periodo)

# 2. Variable binaria resistente/no resistente, a partir de drtype
df_resistoma_resumido$es_resistente <- df_resistoma_resumido$drtype != "Susceptible"
table(df_resistoma_resumido$es_resistente)   # comprobar que se creo bien

# 3. Test de tendencia de Cochran-Armitage -- ¿aumenta la resistencia con el tiempo?
# install.packages("DescTools")  # solo la primera vez
library(DescTools)
tabla_tendencia <- table(df_resistoma_resumido$periodo, df_resistoma_resumido$es_resistente)
CochranArmitageTest(tabla_tendencia, alternative = "increasing")

# 4. Regresion logistica -- controlar pais y periodo a la vez
modelo <- glm(es_resistente ~ isolation_country + periodo,
              data = df_resistoma_resumido, family = binomial)
summary(modelo)
exp(coef(modelo))          # odds ratios
exp(confint(modelo))       # IC 95% de los odds ratios

# 5. Correccion por comparaciones multiples (cuando tengas varios p-valores pais a pais)
p_valores <- c(0.03, 0.01, 0.04, 0.20)   # sustituye por tus p-valores reales
p.adjust(p_valores, method = "BH")

# ==============================================================================
# COMPROBACIÓN ESTADÍSTICA: EFECTO DEL AÑO EN LA RESISTENCIA (CORREGIDO)
# ==============================================================================

# 0. Asegurarnos de que la variable es_resistente existe en el dataframe
df_resistoma_resumido <- df_resistoma_resumido %>%
  mutate(es_resistente = as.integer(drtype != "Susceptible"))

# 1. Prueba de Chi-cuadrado exacta (Monte Carlo)
tabla_temporal <- table(df_resistoma_resumido$collection_year, df_resistoma_resumido$drtype)
test_chi_temporal <- chisq.test(tabla_temporal, simulate.p.value = TRUE, B = 10000)

cat("\n--- RESULTADO CHI-CUADRADO (MONTE CARLO) ---\n")
print(test_chi_temporal)

# 2. Regresión Logística (Binomial)
modelo_logistico <- glm(es_resistente ~ collection_year, 
                        data = df_resistoma_resumido, 
                        family = "binomial")

# ==============================================================================
# SCRIPT: TENDENCIA TEMPORAL DE LA RESISTENCIA, POR PAIS
# ==============================================================================

library(tidyverse)

setwd("C:/Users/migue/Desktop/TFM_R")

df_resistoma_resumido <- read_csv("matriz_resistoma_resumida_muestras_TFM.csv", show_col_types = FALSE) %>%
  mutate(resist_bin = as.integer(drtype != "Susceptible"))

paises <- c("USA", "China", "India", "South Africa", "United Kingdom")

# ------------------------------------------------------------------------
# Regresion logistica independiente por pais: resist_bin ~ year_c (continuo)
# ------------------------------------------------------------------------
resultados_pais <- map_dfr(paises, function(pais) {
  sub <- df_resistoma_resumido %>%
    filter(isolation_country == pais) %>%
    mutate(year_c = collection_year - mean(collection_year))
  
  m <- glm(resist_bin ~ year_c, data = sub, family = binomial)
  
  tibble(
    pais    = pais,
    OR      = exp(coef(m)["year_c"]),
    IC_2.5  = exp(confint(m)["year_c", 1]),
    IC_97.5 = exp(confint(m)["year_c", 2]),
    p       = summary(m)$coefficients["year_c", 4]
  )
})

print(resultados_pais, n = Inf)
write_csv(resultados_pais, "Tendencia_temporal_por_pais.csv")


cat("\n--- RESULTADO REGRESIÓN LOGÍSTICA ---\n")
print(summary(modelo_logistico))

# 3. Calcular el Odds Ratio (OR) y sus Intervalos de Confianza (IC al 95%)
cat("\n--- ODDS RATIO E INTERVALOS DE CONFIANZA ---\n")
resultados_or <- exp(cbind(OR = coef(modelo_logistico), confint(modelo_logistico)))
print(resultados_or)