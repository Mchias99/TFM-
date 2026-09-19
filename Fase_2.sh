#!/bin/bash

# ==========================================================
# PIPELINE AUTORRESCABLE Y ROBUSTO - M. tuberculosis (TFM)
# ==========================================================

WORKDIR="/mnt/c/Users/migue/Desktop/TFM_R"
cd "$WORKDIR" || exit 1

mkdir -p raw_fastq fasta_files tbprofiler_results/results

LOGFILE="pipeline.log"
FAILED_LOG="muestras_fallidas.txt"
exec > >(tee -a "$LOGFILE") 2>&1

echo "=========================================================="
echo "   PIPELINE AUTOMATIZADO DE RESISTOMA CON RESCATE AUTOMÁTICO"
echo "   Iniciado el: $(date)"
echo "=========================================================="

# --------------------------------------------------------
# FASE 1: DESCARGA DE ARCHIVOS SRA
# --------------------------------------------------------
if [ -f "solo_pendientes_sra.txt" ]; then
    sed -i 's/\r$//' solo_pendientes_sra.txt
    echo -e "\n[FASE 1/4] Descargando paquetes SRA..."
    
    while read -r run; do
        [ -z "$run" ] && continue

        if [ -f "tbprofiler_results/results/${run}.results.json" ]; then
            echo "--> [$(date +'%H:%M:%S')] ${run} ya procesada. Omitiendo descarga."
            continue
        fi

        if [ -d "raw_fastq/${run}" ] || [ -f "raw_fastq/${run}/${run}.sra" ]; then
            echo "--> [$(date +'%H:%M:%S')] Archivo SRA para ${run} ya existe en local."
            continue
        fi

        echo "--------------------------------------------------"
        echo "[$(date +'%H:%M:%S')] Descargando paquete SRA para: ${run}"
        prefetch "${run}" -O raw_fastq/ --max-size 50G >/dev/null 2>&1

    done < solo_pendientes_sra.txt
fi

# --------------------------------------------------------
# FASE 2: EXTRACCIÓN + QC + PROFILING + AUTO-RESCATE
# --------------------------------------------------------
if [ -f "solo_pendientes_sra.txt" ]; then
    echo -e "\n=========================================================="
    echo "   [FASE 2/4] Procesamiento y Rescate Automático de Lecturas"
    echo "=========================================================="

    while read -r run; do
        [ -z "$run" ] && continue

        if [ -f "tbprofiler_results/results/${run}.results.json" ]; then
            continue
        fi

        echo "--------------------------------------------------"
        echo "[$(date +'%H:%M:%S')] Procesando muestra: ${run}"

        # 1. Extracción desde SRA usando --split-3 para separar huérfanas
        if [ -f "raw_fastq/${run}/${run}.sra" ] || [ -d "raw_fastq/${run}" ]; then
            echo "[1/3] Extrayendo FASTQ..."
            timeout 15m fasterq-dump --split-3 --skip-technical --outdir raw_fastq --temp /dev/shm "raw_fastq/${run}/${run}.sra" -e 4 >/dev/null 2>&1
            rm -rf "raw_fastq/${run}"
        fi

        SUCCESS=0

        # CASO A: Intento Paired-End estándar
        if [ -f "raw_fastq/${run}_1.fastq" ] && [ -f "raw_fastq/${run}_2.fastq" ]; then
            echo "[2/3] Filtrando lecturas Paired-End con fastp..."
            
            if fastp -i "raw_fastq/${run}_1.fastq" -I "raw_fastq/${run}_2.fastq" \
                      -o "raw_fastq/${run}_1.clean.fastq.gz" -O "raw_fastq/${run}_2.clean.fastq.gz" \
                      --qualified_quality_phred 20 --length_required 50 \
                      --cut_front --cut_tail --thread 8 --html /dev/null --json /dev/null >/dev/null 2>&1; then

                if [ -s "raw_fastq/${run}_1.clean.fastq.gz" ] && [ -s "raw_fastq/${run}_2.clean.fastq.gz" ]; then
                    echo "[3/3] Ejecutando TB-Profiler (Paired-End)..."
                    if tb-profiler profile -1 "raw_fastq/${run}_1.clean.fastq.gz" \
                                           -2 "raw_fastq/${run}_2.clean.fastq.gz" \
                                           -p "${run}" -t 8 --dir tbprofiler_results >/dev/null 2>&1; then
                        SUCCESS=1
                    fi
                fi
            fi

            # MÓDULO DE RESCATE: Si fastp/tb-profiler fallaron por desbalance de pares
            if [ $SUCCESS -eq 0 ]; then
                echo "[RESCATE] Desbalance detectado en ${run}. Ejecutando modo de rescate Single-End..."
                
                shopt -s nullglob
                merge_files=(raw_fastq/"${run}"*.fastq)
                shopt -u nullglob
                if [ "${#merge_files[@]}" -gt 0 ]; then
                    cat "${merge_files[@]}" > "raw_fastq/${run}_merged.fastq"
                fi

                if fastp -i "raw_fastq/${run}_merged.fastq" -o "raw_fastq/${run}_rescued.clean.fastq.gz" \
                          --qualified_quality_phred 20 --length_required 50 \
                          --cut_front --cut_tail --thread 8 --html /dev/null --json /dev/null >/dev/null 2>&1; then

                    if [ -s "raw_fastq/${run}_rescued.clean.fastq.gz" ]; then
                        echo "[RESCATE 3/3] Ejecutando TB-Profiler sobre muestra rescatada..."
                        if tb-profiler profile -1 "raw_fastq/${run}_rescued.clean.fastq.gz" \
                                               -p "${run}" -t 8 --dir tbprofiler_results >/dev/null 2>&1; then
                            echo "[ÉXITO] Muestra ${run} rescatada correctamente."
                            SUCCESS=1
                        fi
                    fi
                fi
            fi

        # CASO B: Single-End desde origen
        elif [ -f "raw_fastq/${run}.fastq" ]; then
            echo "[2/3] Filtrando lecturas Single-End con fastp..."

            if fastp -i "raw_fastq/${run}.fastq" -o "raw_fastq/${run}.clean.fastq.gz" \
                      --qualified_quality_phred 20 --length_required 50 \
                      --cut_front --cut_tail --thread 8 --html /dev/null --json /dev/null >/dev/null 2>&1; then

                if [ -s "raw_fastq/${run}.clean.fastq.gz" ]; then
                    echo "[3/3] Ejecutando TB-Profiler (Single-End)..."
                    if tb-profiler profile -1 "raw_fastq/${run}.clean.fastq.gz" \
                                           -p "${run}" -t 8 --dir tbprofiler_results >/dev/null 2>&1; then
                        SUCCESS=1
                    fi
                fi
            fi
        fi

        # Registrar en log de fallidas si tras el rescate no se procesó
        if [ $SUCCESS -eq 0 ]; then
            echo "[ERROR DEFINITIVO] La muestra ${run} no pudo ser procesada ni rescatada." | tee -a "$FAILED_LOG"
        fi

        # Limpieza de archivos temporales
        rm -rf "raw_fastq/${run}"*.fastq "raw_fastq/${run}"*.fastq.gz

        if [ "$SUCCESS" -eq 1 ]; then
            rm -f "tbprofiler_results/bam/${run}"*.bam "tbprofiler_results/bam/${run}"*.bam.bai
            rm -f "tbprofiler_results/vcf/${run}"*.vcf.gz "tbprofiler_results/vcf/${run}"*.vcf.gz.tbi
        fi

    done < solo_pendientes_sra.txt
fi

# --------------------------------------------------------
# FASE 3: ENSAMBLADOS (FASTA) Y PROCESAMIENTO NCBI
# --------------------------------------------------------
if [ -f "solo_pendientes_gca.txt" ]; then
    sed -i 's/\r$//' solo_pendientes_gca.txt
    echo -e "\n[FASE 3/4] Procesando genomas ensamblados NCBI (FASTA)..."

    while read -r gca; do
        [ -z "$gca" ] && continue

        if [ -f "tbprofiler_results/results/${gca}.results.json" ]; then
            echo "--> [$(date +'%H:%M:%S')] Ensamblado ${gca} ya procesado. Omitiendo."
            continue
        fi

        echo "--------------------------------------------------"
        echo "[$(date +'%H:%M:%S')] Procesando Ensamblado: ${gca}"

        ZIP_FILE="fasta_files/${gca}.zip"
        TEMP_DIR="fasta_files/temp_${gca}"

        timeout 10m datasets download genome accession "${gca}" --filename "${ZIP_FILE}" 2>/dev/null || true

        if [ -f "${ZIP_FILE}" ]; then
            unzip -q -o "${ZIP_FILE}" -d "${TEMP_DIR}" 2>/dev/null || true
            FASTA_FILE=$(find "${TEMP_DIR}" -name "*.fna" -o -name "*.fasta" 2>/dev/null | head -n 1)

            if [ -n "$FASTA_FILE" ] && [ -s "$FASTA_FILE" ]; then
                tb-profiler profile -f "${FASTA_FILE}" -p "${gca}" --dir tbprofiler_results -t 8 || echo "[ERROR] Falló TB-Profiler FASTA en ${gca}"
            fi

            rm -rf "${TEMP_DIR}" "${ZIP_FILE}"
        fi

    done < solo_pendientes_gca.txt
fi

# --------------------------------------------------------
# FASE 4: MÓDULO DE RECUPERACIÓN PARA MUESTRAS FALTANTES
# --------------------------------------------------------
RESCUE_LOG="rescate_resultados.log"

echo -e "\n=========================================================="
echo "   INICIANDO MÓDULO DE RECUPERACIÓN PARA MUESTRAS FALTANTES"
echo "   Iniciado el: $(date)"
echo "=========================================================="

MISSING_SRA="faltantes_sra.txt"
MISSING_GCA="faltantes_gca.txt"
> "$MISSING_SRA"
> "$MISSING_GCA"

# 1. Filtrar exactamente qué muestras SRA faltan por procesar
if [ -f "solo_pendientes_sra.txt" ]; then
    while read -r run; do
        [ -z "$run" ] && continue
        if [ ! -f "tbprofiler_results/results/${run}.results.json" ]; then
            echo "$run" >> "$MISSING_SRA"
        fi
    done < solo_pendientes_sra.txt
fi

# 2. Filtrar exactamente qué ensamblados GCA faltan por procesar
if [ -f "solo_pendientes_gca.txt" ]; then
    while read -r gca; do
        [ -z "$gca" ] && continue
        if [ ! -f "tbprofiler_results/results/${gca}.results.json" ]; then
            echo "$gca" >> "$MISSING_GCA"
        fi
    done < solo_pendientes_gca.txt
fi

NUM_SRA=$(wc -l < "$MISSING_SRA")
NUM_GCA=$(wc -l < "$MISSING_GCA")

echo "Muestras SRA pendientes de rescate: ${NUM_SRA}"
echo "Ensamblados GCA pendientes de rescate: ${NUM_GCA}"
echo "--------------------------------------------------"

# FASE 4A: RESCATE SRA VIA FASTQ-DUMP DIRECTO / ENA VIA HTTP
if [ "$NUM_SRA" -gt 0 ]; then
    while read -r run; do
        [ -z "$run" ] && continue
        echo "[$(date +'%H:%M:%S')] Intentando rescate intensivo para: ${run}"

        rm -rf "raw_fastq/${run}"*

        echo "  -> [1/3] Descargando mediante fastq-dump directo..."
        timeout 20m fastq-dump --split-3 --gzip --outdir raw_fastq "${run}" >/dev/null 2>&1

        if [ ! -f "raw_fastq/${run}_1.fastq.gz" ] && [ ! -f "raw_fastq/${run}.fastq.gz" ]; then
            echo "  -> [2/3] Reintentando vía ENA (FTP Directo)..."
            PREFIX=$(echo "$run" | cut -c 1-6)
            
            curl -sf -L "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${run}/${run}_1.fastq.gz" -o "raw_fastq/${run}_1.fastq.gz" 2>/dev/null
            curl -sf -L "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${run}/${run}_2.fastq.gz" -o "raw_fastq/${run}_2.fastq.gz" 2>/dev/null

            for f in "raw_fastq/${run}_1.fastq.gz" "raw_fastq/${run}_2.fastq.gz"; do
                if [ -s "$f" ] && ! gzip -t "$f" 2>/dev/null; then
                    echo "  [AVISO] Descarga corrupta/truncada detectada en $(basename "$f"), descartando."
                    rm -f "$f"
                elif [ ! -s "$f" ]; then
                    rm -f "$f"
                fi
            done
        fi

        SUCCESS=0

        if [ -f "raw_fastq/${run}_1.fastq.gz" ] && [ -f "raw_fastq/${run}_2.fastq.gz" ]; then
            echo "  -> Descarga exitosa. Filtrando con fastp..."
            if fastp -i "raw_fastq/${run}_1.fastq.gz" -I "raw_fastq/${run}_2.fastq.gz" \
                      -o "raw_fastq/${run}_1.clean.fastq.gz" -O "raw_fastq/${run}_2.clean.fastq.gz" \
                      --qualified_quality_phred 20 --length_required 50 \
                      --cut_front --cut_tail --thread 8 --html /dev/null --json /dev/null >/dev/null 2>&1; then

                tb-profiler profile -1 "raw_fastq/${run}_1.clean.fastq.gz" \
                                    -2 "raw_fastq/${run}_2.clean.fastq.gz" \
                                    -p "${run}" -t 8 --dir tbprofiler_results >/dev/null 2>&1 && SUCCESS=1
            fi
        
        elif [ -f "raw_fastq/${run}.fastq.gz" ] || [ -f "raw_fastq/${run}_1.fastq.gz" ]; then
            SE_FILE=$(ls raw_fastq/${run}*.fastq.gz 2>/dev/null | head -n 1)
            echo "  -> Descarga exitosa (Single-End). Filtrando con fastp..."
            if fastp -i "${SE_FILE}" -o "raw_fastq/${run}_rescued.clean.fastq.gz" \
                      --qualified_quality_phred 20 --length_required 50 \
                      --cut_front --cut_tail --thread 8 --html /dev/null --json /dev/null >/dev/null 2>&1; then

                tb-profiler profile -1 "raw_fastq/${run}_rescued.clean.fastq.gz" \
                                    -p "${run}" -t 8 --dir tbprofiler_results >/dev/null 2>&1 && SUCCESS=1
            fi
        fi

        if [ $SUCCESS -eq 1 ]; then
            echo "  [¡ÉXITO!] Muestra ${run} rescatada y procesada correctamente."
            rm -f "tbprofiler_results/bam/${run}"*.bam "tbprofiler_results/bam/${run}"*.bam.bai
            rm -f "tbprofiler_results/vcf/${run}"*.vcf.gz "tbprofiler_results/vcf/${run}"*.vcf.gz.tbi
        else
            echo "  [FALLO PERMANENTE] No fue posible obtener datos públicos válidos para ${run}."
        fi

        rm -rf "raw_fastq/${run}"*

    done < "$MISSING_SRA"
fi

# FASE 4B: RESCATE GCA (FASTA) VIA REINTENTO DIRECTO
if [ "$NUM_GCA" -gt 0 ]; then
    while read -r gca; do
        [ -z "$gca" ] && continue
        echo "[$(date +'%H:%M:%S')] Reintentando descarga FASTA para: ${gca}"

        ZIP_FILE="fasta_files/${gca}.zip"
        TEMP_DIR="fasta_files/temp_${gca}"

        datasets download genome accession "${gca}" --filename "${ZIP_FILE}" >/dev/null 2>&1

        if [ -f "${ZIP_FILE}" ]; then
            unzip -q -o "${ZIP_FILE}" -d "${TEMP_DIR}" 2>/dev/null || true
            FASTA_FILE=$(find "${TEMP_DIR}" -name "*.fna" -o -name "*.fasta" 2>/dev/null | head -n 1)

            if [ -n "$FASTA_FILE" ] && [ -s "$FASTA_FILE" ]; then
                tb-profiler profile -f "${FASTA_FILE}" -p "${gca}" --dir tbprofiler_results -t 8 >/dev/null 2>&1
                echo "  [¡ÉXITO!] Ensamblado ${gca} procesado correctamente."
            fi
            rm -rf "${TEMP_DIR}" "${ZIP_FILE}"
        fi
    done < "$MISSING_GCA"
fi

awk '/INICIANDO MÓDULO DE RECUPERACIÓN PARA MUESTRAS FALTANTES/{found=1} found' "$LOGFILE" > "$RESCUE_LOG"

# --------------------------------------------------------
# FASE 5: CONSOLIDACIÓN GENERAL FINAL
# --------------------------------------------------------
echo -e "\n=========================================================="
echo "   CONSOLIDANDO MATRIZ GENERAL DE RESULTADOS AL FINAL"
echo "=========================================================="

# Buscamos explícitamente archivos terminados en .results.json
TOTAL_JSON=$(ls -1 tbprofiler_results/results/*.results.json 2>/dev/null | wc -l)

if [ "$TOTAL_JSON" -gt 0 ]; then
    echo "Encontrados $TOTAL_JSON archivos .results.json para consolidar..."
    
    # 1. Entramos al directorio raíz de tb-profiler para evitar bugs de rutas
    cd tbprofiler_results || exit 1
    
    # 2. Ejecutamos collate sin usar --dir. Por defecto buscará en ./results/
    if tb-profiler collate --prefix tbprofiler --full; then
        echo "[ÉXITO] Consolidación finalizada correctamente."
    else
        echo "[ERROR] tb-profiler collate falló internamente."
    fi
    
    # 3. Volvemos al directorio de trabajo principal
    cd ..
else
    echo "[AVISO] No se encontraron archivos .results.json en tbprofiler_results/results/"
fi

# --------------------------------------------------------
# FASE 6: REVISIÓN DE COBERTURA Y MARCADO DE MUESTRAS (FLAGGING)
# --------------------------------------------------------
echo -e "\n=========================================================="
echo "   REVISIÓN DE COBERTURA Y MARCADO DE MUESTRAS (FLAGGING)"
echo "=========================================================="

MAIN_SUMMARY="tbprofiler_results/tbprofiler.txt"
LOW_COV_REPORT="muestras_baja_cobertura_flagged.txt"

if [ -f "$MAIN_SUMMARY" ]; then
    echo "Analizando métricas de cobertura en $MAIN_SUMMARY..."

    awk -F'\t' '
    NR==1 {
        for (i=1; i<=NF; i++) {
            h=$i; gsub(/^[ \t]+|[ \t]+$/, "", h); h=tolower(h)
            if (h ~ /pct.*mapped|percent.*mapped|mapped.*pct/) col_mapped=i
            else if (!col_mapped && h ~ /mapped/) col_mapped_fallback=i
            if (h ~ /median.*(coverage|depth)|mean.*(coverage|depth)/) col_depth=i
            else if (!col_depth && h ~ /coverage|depth/) col_depth_fallback=i
        }
        if (!col_mapped && col_mapped_fallback) {
            print "AVISO: no se encontro columna de PORCENTAJE de mapeado (pct_reads_mapped); usando " col_mapped_fallback " como aproximacion." > "/dev/stderr"
            col_mapped = col_mapped_fallback
        }
        if (!col_depth && col_depth_fallback) col_depth = col_depth_fallback
        if (!col_mapped || !col_depth) {
            print "AVISO: no se localizaron las columnas de %mapeado/cobertura por nombre; revisa la cabecera de " FILENAME > "/dev/stderr"
        }
        next
    }
    col_mapped && col_depth {
        gsub(/^[ \t]+|[ \t]+$/, "", $1);
        mapped = $col_mapped + 0;
        depth  = $col_depth + 0;
        if (mapped > 100.0) {
            printf "Muestra: %-15s | COLUMNA_MAPEADO_INVALIDA (valor=%.2f) | Cobertura: %6.2fx\n", $1, mapped, depth
            next
        }
        if (mapped < 50.0 || depth < 10.0) {
            printf "Muestra: %-15s | %% Mapeado: %6.2f%% | Cobertura: %6.2fx | FLAG: BAJA_COBERTURA\n", $1, mapped, depth
        }
    }' "$MAIN_SUMMARY" > "$LOW_COV_REPORT"

    TOTAL_LOW=$(wc -l < "$LOW_COV_REPORT")

    if [ "$TOTAL_LOW" -gt 0 ]; then
        echo -e "\n[ALERTA QC] Se detectaron ${TOTAL_LOW} muestras con baja cobertura/calidad:"
        echo "----------------------------------------------------------------------"
        cat "$LOW_COV_REPORT"
        echo "----------------------------------------------------------------------"
        echo "El listado completo de alertas se guardó en: ${LOW_COV_REPORT}"
    else
        echo "[QC OK] Todas las muestras procesadas superan los umbrales de calidad (>=50% mapped, >=10x depth)."
    fi
else
    echo "[ERROR QC] No se encontró el resumen consolidado $MAIN_SUMMARY."
fi

echo -e "\n=========================================================="
echo " PIPELINE AUTOMATIZADO COMPLETO FINALIZADO CON ÉXITO"
echo " Finalizado el: $(date)"
if [ -f "$FAILED_LOG" ]; then
    echo " Muestras no rescatables: $(wc -l < "$FAILED_LOG") (ver $FAILED_LOG)"
fi
echo "=========================================================="