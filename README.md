# GEVIS - Differential Gene Expression Analysis with Visual Analytics

Differential gene expression (DGE) analysis is one of the most used techniques for RNA-seq data analysis, and it is applied in various medical and biological contexts, including biomarkers for diagnosis and prognosis and evaluation of the effectiveness of specific treatments. The conduction of a DGE analysis typically involves navigating a complex, multi-step pipeline, which usually requires proficiency in programming languages like R. This presents a barrier to researchers like biologists and clinicians, who may have limited or no coding skills, and adds additional overhead even for experienced bioinformaticians. 

To overcome these challenges, we propose **GEVIS**, visual analytics prototype that simplifies DGE analysis, enabling users to perform the analyses without coding expertise.

## Citation
If you use GEVIS in your research, please cite the following paper:


> Differential Gene Expression Analysis with Visual Analytics.
_Francesco Fortunato, Cristian Santaroni, Graziano Blasilli, Giulia Fiscon, Simone Lenti, Giuseppe Santucci._ In EuroVis 2025 - Posters. The Eurographics Association, 2025. DOI: [10.2312/evp.20251126](https://doi.org/10.2312/evp.20251126).
```
@inproceedings{10.2312:evp.20251126,
  booktitle = {EuroVis 2025 - Posters},
  editor = {Diehl, Alexandra and Kucher, Kostiantyn and Médoc, Nicolas},
  title = {{Differential Gene Expression Analysis with Visual Analytics}},
  author = {Fortunato, Francesco and Santaroni, Cristian and Blasilli, Graziano and Fiscon, Giulia and Lenti, Simone and Santucci, Giuseppe},
  year = {2025},
  publisher = {The Eurographics Association},
  ISBN = {978-3-03868-286-8},
  DOI = {10.2312/evp.20251126}
}
```

<img src="pictures\GEVIS UI 2.png" alt="Dashboard Image" style="max-width: 100%; height: auto;">


<img src="pictures\GEVIS UI 1.png" alt="Dashboard Image" style="max-width: 100%; height: auto;">


# Table of Contents

- [Features](#features)
- [Running Locally](#running-locally)
- [Usage](#usage)
- [Input File Format](#input-file-format)
  - [Matrix File](#matrix-file)
  - [Metadata File](#metadata-file)
- [License](#license)


---

## Features

- **Scatterplot and Volcano Plot Visualization:** 
  - Easily switch between scatterplot (comparing case expression vs. normal expression) and volcano plot for intuitive visual exploration of differentially expressed genes.
  - Scatterplot circles represent genes, color-coded based on log fold change (LogFC). Clicking a gene highlights it and adds its axis to the scatterplot.
  - Volcano plot visualizes gene significance concerning LogFC and p-value, helping identify the most significant genes at a glance.

- **Custom Dataset Input:** 
  - Users can input raw gene expression data along with metadata, offering flexibility in the types of data analyzed. This feature allows for the integration of user-specific datasets for personalized exploration.

- **Multiple Example Datasets:** 
  - Choose from multiple preloaded example datasets for quick analysis, enabling comparative studies across different experiments and conditions.

- **Quantile normalization:** 
  - Perform quantile normalization.

- **Statistical Testing Options:** 
  - Perform differential gene expression analysis using three statistical methods: t-test, Limma, DESeq2, and Wilcoxon. These options offer flexibility depending on the nature of the dataset and the analysis goals.
  - Results include p-values and adjusted p-values (False Discovery Rate - FDR) for identifying statistically significant genes.

- **Enrichment Analysis:** 
  - Conduct enrichment analysis on differentially expressed genes, helping identify pathways, functions, or gene ontologies which are overrepresented in the dataset.
  
- **Survival Analysis:** 
  - Integrate survival data and perform survival analysis based on gene expression levels. This feature is crucial for linking gene expression patterns to patient outcomes.

- **Gene Export:** 
  - Users can export the list of genes identified as significant during the analysis for downstream processes or reporting.

- **Parallel Coordinates Plot with Custom Color Encoding:** 
  - Visualize gene expression profiles across samples using an interactive parallel coordinates plot. The user can select metadata axes, reorder axes, and highlight specific samples.
  - The user can also customize the color encoding of the plot for deeper insights into patterns across different sample groups.

- **Principal Component Analysis (PCA):** 
  - Perform PCA to explore clustering and variability among samples. This feature helps in identifying patterns of gene expression that contribute to sample separation.

- **Heatmap of PCA Contributions:** 
  - Visualize the contributions of genes to the principal components through a heatmap, making it easier to understand which genes drive the variability among samples.

- **Boxplot Visualization:** 
  - Display the distribution of gene expression values across different sample groups using boxplots. This feature allows users to assess variability, outliers, and trends within the data.

---

## Running Locally

1. Open your terminal.  
2. Clone the repository to your local machine by running:  
    ```sh
    git clone https://github.com/francesco-fortunato/GEVIS.git
    ```  
3. Navigate to the "GEVIS" directory using the command line:  
    ```sh
    cd GEVIS
    ```  
4. Install Git LFS (if it is not already installed):  
    ```sh
    git lfs install
    ```  
5. Pull large files managed by Git LFS:  
    ```sh
    git lfs pull
    ```  
6. Build and run the application using Docker Compose:  
    ```sh
    docker-compose pull
    docker-compose up --build
    ```  
7. Open your web browser and go to [localhost:11765](http://localhost:11765) to access the GEVIS dashboard.  
   Use the interactive sliders and graphs to explore gene expression data and identify differentially expressed genes.  

P.S. Examples one and two are too heavy for the heatmap; use example three for a smooth experience. Enjoy :)
---

## Usage

1. **Data Upload:**
   - Upload raw gene expression data with associated metadata or select from available example datasets for analysis.
   
2. **Select Analysis Type:**
   - Choose between t-test, Limma, DESeq, or Wilcoxon for differential expression analysis.

3. **Interactive Visualization:**
   - Explore the results visually through scatterplots, volcano plots, parallel coordinate plots, and PCA. 
   - Customize axes, metadata, and color encoding for tailored visualizations.
   
4. **Perform Enrichment & Survival Analysis:**
   - Utilize built-in tools for enrichment analysis and survival analysis based on the results of your differential expression studies.
   
5. **Export Results:**
   - Export genes and analysis results for further exploration or reporting.

---

## Input File Format

The application requires two input files in CSV format: a **matrix file** and a **metadata file**. Below is the detailed structure for each file:

### Matrix File

The matrix file should be structured as follows:

- The first column must contain the gene identifiers, and the header for this column can have any name.
- The subsequent columns should contain expression values corresponding to each sample, with the column headers representing the sample identifiers.

**Example:**

#### Gene Expression Matrix

| ID_REF         | Sample 1           | Sample 2           | Sample 3           | ... |
|----------------|----------------------|----------------------|----------------------|-----|
| Gene A | 230     | 300     | 240     | ... |
| Gene B | 518     | 551     | 569     | ... |
| ...    | ...     | ...     | ...     | ... |

### Metadata File

The metadata file should be structured as follows:

- The first column must contain metadata labels, which can have any descriptive name. This column may include various types of information about each sample, such as disease type, treatment group, or other relevant characteristics.
- Each subsequent column should correspond to a specific sample, containing the relevant metadata for that sample.

**Example:**

#### Sample Metadata

| Metadata | Sample 1           | Sample 2           | Sample 3           | ... |
|----------------------|----------------------|----------------------|----------------------|-----|
| Type                 | Intrahepatic cholangiocarcinoma | Intrahepatic cholangiocarcinoma | Normal intrahepatic bile duct | ... |
| Gender               | Male                 | Female               | Male                 | ... |
| ...                  | ...                  | ...                  | ...                  | ... |

Ensure that both files are formatted correctly to conduct proper analysis within the application.
---

## License

This software is provided under a strict proprietary license. All rights to the software, including its code and associated files, are retained by the author. Unauthorized use, reproduction, or redistribution of any part of the software is prohibited. The software is made available solely for educational or review purposes. Any commercial or derivative use requires explicit written permission from the author.

For licensing inquiries, please contact: [francesco.fortunato1999@gmail.com](mailto:francesco.fortunato1999@gmail.com).
