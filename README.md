# Single-Cell RNA-seq Analysis Pipeline

A **Nextflow-based workflow** for processing and analysis of **10x Genomics single-cell RNA-seq data**.

## Workflow overview

This pipeline performs comprehensive scRNA-seq data analysis from raw data preprocessing through cell type annotation. All tools run in isolated Docker containers, ensuring reproducibility across different computing environments.

---

## Workflow Summary

### 1. Read Quality Assessment
- Quality control using **FastQC** (v0.12.1)
- Generates HTML reports for each sample (R1 and R2)

### 2. Alignment and Quantification
- Aligns reads to reference transcriptome using **Cell Ranger** (v9.0.1)
- Outputs:
  - `web_summary.html` — alignment metrics
  - `filtered_feature_bc_matrix.h5` — gene-barcode matrix

### 3. MultiQC Report
- Aggregates QC metrics from FastQC and Cell Ranger using **MultiQC** (v1.21)

### 4. Scanpy Analysis

**Scanpy** (v1.11.4) performs downstream single-cell analysis:

#### 4.1 Quality Control and Filtering
- Cell filtering: min 500 UMI, <10% mitochondrial genes, MAD outlier detection
- Doublet detection: Scrublet (threshold = 0.25)
- Gene filtering: expressed in ≥3 cells

#### 4.2 Normalization and Feature Selection
- Log-normalization (target_sum = 10,000)
- 4,000 highly variable genes (cell_ranger flavor)

#### 4.3 Batch Correction
- Harmony integration across samples

#### 4.4 Dimensionality Reduction
- PCA: 50 components → 20 used
- UMAP visualization

#### 4.5 Clustering
- Leiden algorithm with multiple resolutions (0.3, 0.5, 0.8, 1.2)
- Parameters: 15 neighbors, 20 PCs

#### 4.6 Cell Type Annotation
Marker genes visualization:
- T cells: Cd3d, Cd3e
- B cells: Cd19, Ms4a1
- NK cells: Ncr1, Nkg7
- Monocytes: Itgam, Ly6c2
- Dendritic cells: Itgax, Cd209a
- Plasma cells: Sdc1, Jchain

#### 4.7 Outputs
- `adata_processed.h5ad` — processed AnnData object
- `cell_metadata.csv` — cell metadata
- `cluster_counts.csv` — cell counts per cluster
- `qc_summary.csv` — QC statistics
- Multiple visualization plots (PNG)

---

## Usage

### Quick Start

```bash
nextflow run main.nf
```

**Resume interrupted run:**
```bash
nextflow run main.nf -resume
```

---

## Options (Flags)

### Analysis mode

- `--skip_scanpy`  
  Skip Scanpy analysis module (only run FastQC, Cell Ranger, and MultiQC)

**Skip Scanpy analysis:**
```bash
nextflow run main.nf --skip_scanpy
```

---

## Requirements

- **[Nextflow](https://www.nextflow.io/docs/latest/install.html)**

- **[Docker](https://docs.docker.com/engine/install/)**

- **Reference transcriptome**
  - Mouse reference: `refdata-gex-mm10-2020-A` from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
  - Pre-built STAR index and gene annotations included

- **FASTQ files**
  - Format: `{sample}_S1_L001_R{1,2}_001.fastq.gz`
  - Standard 10x Genomics naming convention
  - Both R1 (Read 1) and R2 (Read 2) required for each sample

---

## Running the Pipeline

1. **Install Nextflow and Docker**

2. **Clone the repository:**
   ```bash
   git clone https://github.com/genespace-workflows/sc_analysis.git
   ```

3. **Navigate to the pipeline directory:**
   ```bash
   cd sc_analysis
   ```

4. **Download reference transcriptome:**
   ```bash
   wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
   tar -xzf refdata-gex-mm10-2020-A.tar.gz
   ```

5. **Edit `nextflow.config`** to match your data paths


6. **Run the pipeline:**
   ```bash
   nextflow run main.nf
   ```

---


## License

*This project is licensed under the MIT License.*
