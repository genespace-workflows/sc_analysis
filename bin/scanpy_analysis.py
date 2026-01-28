#!/usr/bin/env python3
"""
Скрипт для анализа single-cell RNA-seq данных с помощью Scanpy
"""

import argparse
import os
import glob
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import median_abs_deviation
from pathlib import Path
import anndata as ad

# Настройки Scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(6, 6))

def parse_args():
    """Парсинг аргументов командной строки"""
    parser = argparse.ArgumentParser(description='Анализ single-cell данных в Scanpy')

    # Входные/выходные данные
    parser.add_argument('--input_dir', type=str, required=True,
                        help='Директория с .h5 файлами от Cell Ranger')
    parser.add_argument('--output_dir', type=str, default='scanpy_results',
                        help='Директория для сохранения результатов')

    # Параметры фильтрации
    parser.add_argument('--min_counts', type=int, default=500,
                        help='Минимальное число UMIs на клетку')
    parser.add_argument('--min_cells', type=int, default=3,
                        help='Минимальное число клеток для гена')
    parser.add_argument('--mito_threshold', type=float, default=10.0,
                        help='Порог митохондриальных генов (%)')
    parser.add_argument('--doublet_threshold', type=float, default=0.25,
                        help='Порог для Scrublet doublet score')
    parser.add_argument('--nmads', type=int, default=3,
                        help='MAD множитель для детекции выбросов')

    # Параметры анализа
    parser.add_argument('--n_top_genes', type=int, default=4000,
                        help='Количество высоковариабельных генов')
    parser.add_argument('--n_pcs', type=int, default=20,
                        help='Количество главных компонент для neighbors')
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help='Количество ближайших соседей')
    parser.add_argument('--leiden_resolution', type=float, default=1.2,
                        help='Резолюция для финальной кластеризации Leiden')

    return parser.parse_args()

def plot_qc_distributions(adata, sample_id, output_dir):
    """Построение графиков распределения QC метрик"""
    fig, axs = plt.subplots(1, 4, figsize=(16, 4))

    # 1. Распределение UMI (log1p UMI per cell)
    axs[0].hist(adata.obs['log1p_total_counts'], bins=50, color='blue', edgecolor='black', linewidth=0.7)
    axs[0].set_title('UMI Distribution')
    axs[0].set_xlabel('log1p(UMI per cell)')
    axs[0].set_ylabel('Number of cells')

    # 2. Распределение экспрессированных генов
    axs[1].hist(adata.obs['n_genes_by_counts'], bins=50, color='green', edgecolor='black', linewidth=0.7)
    axs[1].set_title('Distribution of Expressed Genes')
    axs[1].set_xlabel('Number of expressed genes per cell')
    axs[1].set_ylabel('Number of cells')

    # 3. Распределение генов по количеству клеток
    gene_counts = np.sum(adata.X > 0, axis=0).A1 if hasattr(adata.X, 'A1') else np.sum(adata.X > 0, axis=0)
    axs[2].hist(gene_counts, bins=50, log=True, color="purple", edgecolor="black", linewidth=0.7) 
    axs[2].set_title("Gene detection distribution") 
    axs[2].set_xlabel('Number of cells expressing a gene') 
    axs[2].set_ylabel('Number of genes')

    # 4. Распределение митохондриальных генов
    axs[3].hist(adata.obs['pct_counts_mt'], bins=50, color='salmon', edgecolor='black', linewidth=0.7)
    axs[3].set_title('Distribution of MT Genes')
    axs[3].set_xlabel('% of mitochondrial genes per cell')
    axs[3].set_ylabel('Number of cells')

    fig.suptitle(f'{sample_id} - QC Distributions')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_id}_qc_distributions.png'), dpi=150, bbox_inches='tight')
    plt.close()

def is_outlier(adata, metric: str, nmads: int):
    """Определение выбросов по MAD (Median Absolute Deviation)"""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def load_and_prepare_data(input_dir):
    """Загрузка всех .h5 файлов с метаданными"""
    h5_files = glob.glob(os.path.join(input_dir, '*.h5'))

    if len(h5_files) == 0:
        raise FileNotFoundError(f"Не найдено .h5 файлов в {input_dir}")

    print(f"Найдено {len(h5_files)} образцов")

    # Определяем sample names из имен файлов
    sample_dirs = []
    for h5_file in sorted(h5_files):
        # Извлекаем имя образца из пути
        sample_name = Path(h5_file).stem
        # Убираем префикс 'filtered_feature_bc_matrix' если есть
        if 'filtered_feature_bc_matrix' in sample_name:
            sample_name = Path(h5_file).parent.parent.parent.name
        sample_dirs.append(sample_name)

    # Создаем метаданные
    metadata = pd.DataFrame({
        'sample': sample_dirs,
        'condition': ['B_1_351' if 'B_1_351' in s else 'B_1_617' if 'B_1_617' in s else 'PBS' if 'PBS' in s else 'WA1' for s in sample_dirs],
        'concentration': ['10' if '10_' in s else '1' if '1_' in s else 'NA' for s in sample_dirs]
    })

    return h5_files, metadata

def quality_control_single_sample(adata, sample_id, args, output_dir):
    """Контроль качества для одного образца"""
    print(f'\nProcessing {sample_id}:')
    print(f"Total number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")

    # Сохраняем исходные значения
    initial_n_cells = adata.n_obs
    initial_n_genes = adata.n_vars

    # Помечаем митохондриальные и рибосомальные гены
    adata.var["mt"] = adata.var_names.str.startswith("mt-")  
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))

    # Рассчитываем метрики качества
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        inplace=True,
        percent_top=None,
        log1p=True
    )

    # Построение графиков до фильтрации
    plot_qc_distributions(adata, sample_id, output_dir)

    # Фильтруем клетки с низкой глубиной прочтения
    sc.pp.filter_cells(adata, min_counts=args.min_counts)

    # Выявление выбросов по MAD
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", args.nmads)
        | is_outlier(adata, "log1p_n_genes_by_counts", args.nmads)
    )
    print(f'Outlier counts: {adata.obs.outlier.sum()} / {adata.n_obs}')

    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", args.nmads) | (adata.obs["pct_counts_mt"] > args.mito_threshold)
    print(f'MT outlier counts: {adata.obs.mt_outlier.sum()} / {adata.n_obs}')

    # Фильтруем найденные выбросы
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    # Фильтруем гены
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Number of genes after gene filter: {adata.n_vars}")

    # Постфильтрационные графики
    plot_qc_distributions(adata, f"{sample_id}_filtered", output_dir)
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    # Детекция дублетов с Scrublet
    try:
        import scrublet as scr
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        adata.obs["doublet_score"] = doublet_scores
        adata.obs["doublet"] = predicted_doublets
        adata = adata[adata.obs["doublet_score"] < args.doublet_threshold].copy()
        print(f"Count of predicted doublets removed: {predicted_doublets.sum()}")

        # График распределения doublet scores
        plt.figure(figsize=(8, 4))
        plt.hist(doublet_scores, bins=50, edgecolor='black')
        plt.xlabel("Doublet score")
        plt.ylabel("Number of cells")
        plt.title(f"{sample_id} - Doublet Score Distribution")
        plt.axvline(args.doublet_threshold, color='red', linestyle='--', label=f'Threshold={args.doublet_threshold}')
        plt.legend()
        plt.savefig(os.path.join(output_dir, f'{sample_id}_doublet_scores.png'), dpi=150, bbox_inches='tight')
        plt.close()

    except ImportError:
        print("WARNING: scrublet not installed, skipping doublet detection")

    # Финальные метрики
    final_n_cells = adata.n_obs
    final_n_genes = adata.n_vars

    cell_retention_pct = (final_n_cells / initial_n_cells * 100) if initial_n_cells > 0 else 0
    gene_retention_pct = (final_n_genes / initial_n_genes * 100) if initial_n_genes > 0 else 0

    summary = {
        "sample": sample_id,
        "initial_cells": initial_n_cells,
        "final_cells": final_n_cells,
        "cell_retention_%": round(cell_retention_pct, 2),
        "initial_genes": initial_n_genes,
        "final_genes": final_n_genes,
        "gene_retention_%": round(gene_retention_pct, 2)
    }

    return adata, summary

def normalize_and_integrate(adata, args, output_dir):
    """Нормализация, отбор генов, PCA и Harmony интеграция"""
    print("\n=== Нормализация и интеграция ===")

    # Histogram до нормализации
    plt.figure(figsize=(8, 4))
    sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
    plt.xlabel("Total counts")
    plt.ylabel("Number of cells")
    plt.title("Total counts before normalization")
    plt.savefig(os.path.join(output_dir, 'counts_before_normalization.png'), dpi=150, bbox_inches='tight')
    plt.close()

    # Сохраняем raw данные
    adata.raw = adata.copy()

    # Нормализация
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Histogram после нормализации
    normalized_sums = np.asarray(adata.X.sum(axis=1)).flatten()
    plt.figure(figsize=(8, 4))
    sns.histplot(normalized_sums, bins=100, kde=False)
    plt.xlabel("Normalized counts")
    plt.ylabel("Number of cells")
    plt.title("Total counts after normalization")
    plt.savefig(os.path.join(output_dir, 'counts_after_normalization.png'), dpi=150, bbox_inches='tight')
    plt.close()

    # Создаём копию для batch-коррекции
    adata_corrected = adata.copy()

    # Выбор высоковариабельных генов
    sc.pp.highly_variable_genes(
        adata_corrected,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=args.n_top_genes,
        flavor="cell_ranger",
        batch_key="sample"
    )

    print(f"Highly variable genes selected: {adata_corrected.var['highly_variable'].sum()}")

    # Масштабирование
    sc.pp.scale(adata_corrected, max_value=10)

    # PCA
    sc.tl.pca(adata_corrected, n_comps=50, use_highly_variable=True)

    # График variance ratio
    sc.pl.pca_variance_ratio(adata_corrected, n_pcs=50, save='_variance_ratio.png')

    # Harmony интеграция
    try:
        import scanpy.external as sce
        sce.pp.harmony_integrate(adata_corrected, "sample", basis='X_pca')
        print("Harmony batch correction completed")
    except ImportError:
        print("WARNING: scanpy.external (harmony) not available, using PCA without batch correction")
        adata_corrected.obsm['X_pca_harmony'] = adata_corrected.obsm['X_pca'].copy()

    return adata_corrected

def clustering_and_visualization(adata, args, output_dir):
    """Кластеризация Leiden и UMAP визуализация"""
    print("\n=== Кластеризация и визуализация ===")

    # Построение kNN графа
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

    # UMAP
    sc.tl.umap(adata)

    # UMAP по sample, condition, concentration
    sc.pl.umap(adata, color=["sample", "condition", "concentration"], wspace=0.5, save='_metadata.png')

    # Leiden кластеризация с разными резолюциями
    resolutions = [0.3, 0.5, 0.8, 1.2]
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
        sc.pl.umap(
            adata,
            color=f'leiden_{res}',
            title=f'Leiden resolution {res}',
            wspace=0.5,
            legend_loc='on data',
            save=f'_leiden_{res}.png'
        )

    # Финальная кластеризация с выбранной резолюцией
    final_res = args.leiden_resolution
    if f'leiden_{final_res}' not in adata.obs.columns:
        sc.tl.leiden(adata, resolution=final_res, key_added=f'leiden_{final_res}')

    print(f"Final clustering (resolution={final_res}): {len(adata.obs[f'leiden_{final_res}'].unique())} clusters")

    return adata

def annotate_and_visualize_markers(adata, args, output_dir):
    """Визуализация маркерных генов"""
    print("\n=== Аннотация маркерными генами ===")

    # Маркерные гены
    genes_of_interest = ["Cd3d", "Cd19", "Ncr1", "Itgam", "Itgax", "Sdc1"]
    genes_present = [g for g in genes_of_interest if g in adata.var_names]

    if len(genes_present) == 0:
        print("WARNING: None of the marker genes found in dataset")
        return adata

    print(f"Marker genes found: {genes_present}")

    # UMAP с кластерами и маркерными генами
    leiden_key = f'leiden_{args.leiden_resolution}'
    sc.pl.umap(
        adata,
        color=[leiden_key] + genes_present,
        cmap='Reds',
        ncols=4,
        frameon=False,
        vmin=0,
        vmax='p99',
        use_raw=True,
        title=['Leiden clusters'] + genes_present,
        legend_loc='on data',
        save='_clusters_and_genes.png'
    )

    # Dotplot
    sc.pl.dotplot(
        adata,
        var_names=genes_present,
        groupby=leiden_key,
        use_raw=True,
        cmap='Reds',
        dot_max=1.0,
        dot_min=0.0,
        standard_scale='var',
        title='Marker gene expression by cluster',
        save='_marker_genes_dotplot.png'
    )

    return adata

def save_results(adata, output_dir, summary_df):
    """Сохранение результатов"""
    print(f"\n=== Сохранение результатов в {output_dir} ===")

    # Сохраняем AnnData объект
    adata.write(os.path.join(output_dir, 'adata_processed.h5ad'))
    print(f"Saved: adata_processed.h5ad ({adata.n_obs} cells, {adata.n_vars} genes)")

    # Таблица с метаданными клеток
    adata.obs.to_csv(os.path.join(output_dir, 'cell_metadata.csv'))
    print("Saved: cell_metadata.csv")

    # Количество клеток по образцам и кластерам
    leiden_cols = [col for col in adata.obs.columns if col.startswith('leiden_')]
    if leiden_cols:
        final_leiden = leiden_cols[-1]
        cluster_counts = adata.obs.groupby(['sample', 'condition', final_leiden]).size().reset_index(name='n_cells')
        cluster_counts.to_csv(os.path.join(output_dir, 'cluster_counts.csv'), index=False)
        print("Saved: cluster_counts.csv")

    # Сводка QC
    summary_df.to_csv(os.path.join(output_dir, 'qc_summary.csv'), index=False)
    print("Saved: qc_summary.csv")

    print("\n" + "="*80)
    print("QC SUMMARY")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80)

def main():
    """Основная функция"""
    args = parse_args()

    print("="*80)
    print("Scanpy Single-Cell Analysis Pipeline")
    print("="*80)
    print(f"Scanpy version: {sc.__version__}")

    # Создаём выходную директорию
    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir

    # 1. Загрузка данных
    print("\n=== Загрузка данных ===")
    h5_files, metadata = load_and_prepare_data(args.input_dir)

    # 2. Контроль качества для каждого образца
    print("\n=== Контроль качества ===")
    filtered_adatas = {}
    summary_data = []

    for h5_file, (_, row) in zip(sorted(h5_files), metadata.iterrows()):
        sample_id = row['sample']

        # Загружаем образец
        adata = sc.read_10x_h5(h5_file)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        # Добавляем метаданные
        adata.obs['sample'] = sample_id
        adata.obs['condition'] = row['condition']
        adata.obs['concentration'] = row['concentration']

        # QC фильтрация
        adata_filtered, summary = quality_control_single_sample(adata, sample_id, args, args.output_dir)

        filtered_adatas[sample_id] = adata_filtered
        summary_data.append(summary)

    summary_df = pd.DataFrame(summary_data)

    # 3. Объединение данных
    print("\n=== Объединение данных ===")
    adata_pooled = ad.concat(
        adatas=filtered_adatas,
        join='outer',
        merge='unique',
        index_unique='-'
    )
    print(f"Pooled data - Cells: {adata_pooled.n_obs}, Genes: {adata_pooled.n_vars}")

    # 4. Нормализация и интеграция
    adata_corrected = normalize_and_integrate(adata_pooled, args, args.output_dir)

    # 5. Кластеризация и UMAP
    adata_corrected = clustering_and_visualization(adata_corrected, args, args.output_dir)

    # 6. Визуализация маркерных генов
    adata_corrected = annotate_and_visualize_markers(adata_corrected, args, args.output_dir)

    # 7. Сохранение результатов
    save_results(adata_corrected, args.output_dir, summary_df)

    print("\n" + "="*80)
    print("АНАЛИЗ ЗАВЕРШЁН УСПЕШНО!")
    print("="*80)

if __name__ == '__main__':
    main()
