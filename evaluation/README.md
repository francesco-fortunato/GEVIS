# Differential Gene Expression Analysis Pipeline

## Workflow Overview
1. **Data Preprocessing**: Script_evaluation.R
2. **Principal Component Analysis**: PCA.R
3. **Functional Enrichment**: Enrichment.R
4. **Survival Analysis**: Survival.R
5. **Interaction Network**: Interaction_network.R

## 1. Data Preprocessing (Script_evaluation.R)
**Purpose**: Differential Gene Expression (DEG) Analysis

### Key Steps:
- Log transformation of data
- IQR filtering
- Fold change calculation
- Statistical analysis (t-test)
- Multiple testing correction
- Visualization:
  - Volcano plot
  - Pie chart of up/down-regulated genes
  - Box plots for top genes
  - Heatmap of gene expression

## 2. Principal Component Analysis (PCA.R)
**Purpose**: Dimensionality Reduction and Sample Clustering

### Key Features:
- Transform high-dimensional gene expression data
- Visualize sample distributions
- Compute variance explained by principal components
- Generate:
  - Score plots
  - Scree/pareto plots
  - Sample clustering visualization

## 3. Functional Enrichment (Enrichment.R)
**Purpose**: Biological Interpretation of Differentially Expressed Genes

### Supported Databases:
- GO Molecular Function
- GO Biological Process
- KEGG Pathways
- DisGeNET

### Key Functions:
- Perform enrichment analysis
- Generate bar and dot plots
- Save detailed enrichment results

## 4. Survival Analysis (Survival.R)
**Purpose**: Correlate Gene Expression with Clinical Outcomes

### Key Features:
- Stratify samples into HIGH/LOW expression groups
- Generate Kaplan-Meier survival curves
- Calculate statistical significance
- Analyze gene expression impact on survival

## 5. Interaction Network (Interaction_network.R)
**Purpose**: Visualize miRNA-Target Gene Interactions

### Key Functions:
- Retrieve target genes for specified miRNAs
- Create interaction network graph
- Visualize network relationships
- Highlight miRNA connections

## Prerequisites
- R (version 3.6+)
- Required packages: multiMiR, igraph, survminer, enrichR, ggplot2, pheatmap, etc.

```
install.packages(c("ggplot2", "pheatmap", "survminer", "enrichR", "ggbiplot"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("multiMiR", "survival"))
```

## Usage
1. Choose input files
2. Update file paths in scripts
3. Run scripts sequentially

## Customization
- Adjust thresholds:
  - IQR cutoff
  - Fold change cutoff
  - P-value thresholds
