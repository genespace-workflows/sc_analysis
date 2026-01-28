#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline для анализа single-cell RNA-seq данных
 * Шаги: FastQC -> Cell Ranger -> MultiQC -> Scanpy
 */

// Вывод параметров запуска
log.info """
=========================================
Single-Cell RNA-seq Analysis Pipeline
=========================================
Входная директория с FASTQ : ${params.fastq_dir}
Референсный транскриптом    : ${params.reference}
Директория результатов      : ${params.outdir}
Запуск Scanpy анализа       : ${params.skip_scanpy ? 'НЕТ' : 'ДА'}
=========================================
"""

// Список образцов (24 образца)
def samples = [
    "B_1_351_10_1", "B_1_351_10_2", "B_1_351_10_3",
    "B_1_351_1_1", "B_1_351_1_2", "B_1_351_1_3",
    "B_1_617_10_1", "B_1_617_10_2", "B_1_617_10_3",
    "B_1_617_1_1", "B_1_617_1_2", "B_1_617_1_3",
    "PBS_1", "PBS_2", "PBS_3", "PBS_4", "PBS_5", "PBS_6",
    "WA1_10_1", "WA1_10_2", "WA1_10_3",
    "WA1_1_1", "WA1_1_2", "WA1_1_3"
]

/*
 * ПРОЦЕСС 1: Контроль качества с FastQC
 */
process FASTQC {
    tag "${sample}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_results

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

/*
 * ПРОЦЕСС 2: Выравнивание с Cell Ranger
 */
process CELLRANGER_COUNT {
    tag "${sample}"
    publishDir "${params.outdir}/cellranger", mode: 'copy'

    container 'litd/docker-cellranger:v9.0.1'

    input:
    tuple val(sample), path(reads)

    output:
    path "${sample}", emit: cellranger_dir
    path "${sample}/outs/web_summary.html", emit: web_summary
    path "${sample}/outs/metrics_summary.csv", emit: metrics
    path "${sample}/outs/filtered_feature_bc_matrix.h5", emit: h5_matrix

    script:
    """
    # Создаём временную директорию для FASTQ
    mkdir -p fastq_input
    ln -s \${PWD}/${reads[0]} fastq_input/
    ln -s \${PWD}/${reads[1]} fastq_input/

    cellranger count \
        --id=${sample} \
        --transcriptome=${params.reference} \
        --fastqs=fastq_input \
        --sample=${sample} \
        --create-bam=true \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()}
    """
}

/*
 * ПРОЦЕСС 3: Агрегация отчётов с MultiQC
 */
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    input:
    path fastqc_files
    path cellranger_summaries

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

/*
 * ПРОЦЕСС 4: Анализ в Scanpy (опциональный)
 */
process SCANPY_ANALYSIS {
    publishDir "${params.outdir}/scanpy", mode: 'copy'

    container 'gcfntnu/scanpy:1.11.4'

    input:
    path h5_files

    output:
    path "scanpy_results", emit: results
    path "scanpy_results/*.h5ad", emit: adata
    path "scanpy_results/*.png", emit: plots, optional: true
    path "scanpy_results/*.csv", emit: tables

    script:
    """
    python ${projectDir}/bin/scanpy_analysis.py \
        --input_dir . \
        --output_dir scanpy_results \
        --min_counts 500 \
        --min_cells 3 \
        --mito_threshold 10.0 \
        --doublet_threshold 0.25 \
        --nmads 3 \
        --n_top_genes 4000 \
        --n_pcs 20 \
        --n_neighbors 15 \
        --leiden_resolution 1.2
    """
}

/*
 * ОСНОВНОЙ WORKFLOW
 */
workflow {

    // Создание канала с парами FASTQ файлов для каждого образца
    Channel
        .fromList(samples)
        .map { sample ->
            def r1 = file("${params.fastq_dir}/${sample}_S1_L001_R1_001.fastq.gz")
            def r2 = file("${params.fastq_dir}/${sample}_S1_L001_R2_001.fastq.gz")

            // Проверка существования файлов
            if (!r1.exists() || !r2.exists()) {
                log.error "ОШИБКА: Не найдены FASTQ файлы для образца ${sample}"
                log.error "Ожидаемые файлы:"
                log.error "  R1: ${r1}"
                log.error "  R2: ${r2}"
                exit 1
            }

            return tuple(sample, [r1, r2])
        }
        .set { fastq_ch }

    // Шаг 1: FastQC
    FASTQC(fastq_ch)

    // Шаг 2: Cell Ranger
    CELLRANGER_COUNT(fastq_ch)

    // Шаг 3: MultiQC (собираем FastQC + Cell Ranger отчёты)
    MULTIQC(
        FASTQC.out.fastqc_results.collect(),
        CELLRANGER_COUNT.out.web_summary.collect()
    )

    // Шаг 4: Scanpy (если нужно пропустить --skip_scanpy)
    if (!params.skip_scanpy) {
        SCANPY_ANALYSIS(
            CELLRANGER_COUNT.out.h5_matrix.collect()
        )
    }
}

workflow.onComplete {
    log.info """
    =========================================
    Pipeline завершён!
    =========================================
    Статус      : ${workflow.success ? 'УСПЕШНО' : 'ОШИБКА'}
    Результаты  : ${params.outdir}
    Время       : ${workflow.duration}
    =========================================
    """
}
