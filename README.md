# RNA-Seq Differential Gene Expression Pipeline (DESeq2)

An automated, reproducible R pipeline for **RNA-Seq Differential Gene Expression (DGE)** analysis using Bioconductor's `DESeq2`. This workflow covers count matrix simulation, quality control, statistical testing, and publication quality visualizations.

## Pipeline Overview
1. **Data Simulation & Preprocessing:** Generates a negative-binomial distributed RNA-Seq count matrix (1,000 genes across 3 Control and 3 Treated replicates) with planted upregulated and downregulated transcripts, followed by low count filtering.
2. **Statistical Modeling (`DESeq2`):** Performs size factor normalization, empirical Bayes dispersion estimation, and Wald testing with Benjamini-Hochberg FDR correction.
3. **Exploratory Quality Control:** Applies Variance Stabilizing Transformation (`VST`) and performs Principal Component Analysis (PCA) to verify sample clustering by biological condition.
4. **Data Visualization:**
   - **PCA Plot:** Evaluates sample variance and batch separation.
   - **Volcano Plot:** Highlights significantly upregulated and downregulated genes ($\text{FDR} < 0.05$, $\vert{}\log_2\text{FC}\vert{} \ge 1$).
   - **Hierarchical Clustering Heatmap:** Displays Z-score normalized expression profiles of the top 30 differentially expressed genes.

## Project Structure
```text
RNA-Seq-DESeq2-Pipeline/
├── rnaseq_deseq2_pipeline.R   # Self-contained R pipeline script
├── results/                   # Automatically generated upon execution
│   ├── differential_expression_results.csv
│   ├── pca_plot.png
│   ├── volcano_plot.png
│   └── top30_degs_heatmap.png
└── README.md                  # Project documentation
```

## Getting Started

### Prerequisites
- R (version 4.1 or higher)
- Required packages (`DESeq2`, `ggplot2`, `pheatmap`, `ggrepel`) are automatically checked and installed by the script on first run.

### Running the Pipeline
Run the script from your terminal or within RStudio:
```bash
Rscript rnaseq_deseq2_pipeline.R
```

## Output Summary
All statistical tables and high-resolution (300 DPI) figures are automatically exported to the `results/` folder:
- `differential_expression_results.csv`: Complete table containing `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj` (FDR).
- `pca_plot.png`: Sample clustering across PC1 and PC2.
- `volcano_plot.png`: Visual summary of fold changes vs. statistical significance.
- `top30_degs_heatmap.png`: Clustered expression heatmap of top candidate biomarkers.
