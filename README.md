# GEVIS - Gene Expression VIS

GEVIS is a comprehensive dashboard for conducting differential gene expression analysis, designed to facilitate interactive visualization and in-depth exploration of gene expression data. This tool supports multiple analysis methods and visualization techniques, providing researchers with a versatile platform for understanding gene expression patterns and identifying key biomarkers across biological conditions.

# Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Download](#download)
- [Usage](#usage)
- [Input File Format](#input-file-format)
  - [Matrix File](#matrix-file)
  - [Metadata File](#metadata-file)
- [Docs](#docs)
- [Contributors](#contributors)
- [License](#license)

## Introduction 

Gene expression analysis is among the most important means of uncovering biological mechanisms of diseases and identifying potential biomarkers. Whereas there are many tools to perform DE analyses, in particular within programming environments like R, there is significant hindrance to overcome given usually limited programming expertise among researchers. To this end, we propose herein GEVIS-a user-friendly, web application that enables access to advanced gene expression analytical techniques. Methods of integrated DE analysis, like t-tests, Limma, and DESeq2 in one GUI, allow users to do their analysis without prior coding. In this way, researchers, biologists, and doctors will be enabled to reduce efforts in interpreting the results instead of writing complex code. It enables the investigation of raw and pre-processed gene expression datasets into a rich set of visualizations, such as histograms, volcano plots, scatterplots, and PCA graphs that enable deep insight. Furthermore, GEVIS contains filtering and normalization functionalities and enrichment and survival analysis functionalities within one comprehensive analytical workflow. In this respect, for overcoming the limitations often observed in stand-alone tools, GEVIS is designed to enable the user to combine multiple analytical approaches and compare results among diverse computational frameworks. This flexibility is crucial when one wants to ensure robust biological interpretations. 

<img src="pictures/GEVIS Dashboard v1.0.png" alt="Dashboard Image" style="max-width: 100%; height: auto;">

---

## Features

- **Scatterplot and Volcano Plot Visualization:** 
  - Easily switch between scatterplot (comparing case expression vs. normal expression) and volcano plot for intuitive visual exploration of differentially expressed genes.
  - Scatterplot circles represent genes, color-coded based on log fold change (LogFC). Clicking a gene highlights it and adds its axis to the scatterplot.
  - Volcano plot visualizes gene significance with respect to LogFC and p-value, helping identify the most significant genes at a glance.

- **Custom Dataset Input:** 
  - Users can input raw gene expression data along with metadata, offering flexibility in the types of data analyzed. This feature allows for the integration of user-specific datasets for personalized exploration.

- **Multiple Example Datasets:** 
  - Choose from multiple preloaded example datasets for quick analysis, enabling comparative studies across different experiments and conditions.

- **Quantile normalization:** 
  - Perform quantile normalization.

- **Statistical Testing Options:** 
  - Perform differential gene expression analysis using three statistical methods: t-test, Limma, and DESeq. These options offer flexibility depending on the nature of the dataset and the analysis goals.
  - Results include p-values and adjusted p-values (False Discovery Rate - FDR) for identifying statistically significant genes.

- **Enrichment Analysis:** 
  - Conduct enrichment analysis on differentially expressed genes, helping identify pathways, functions, or gene ontologies that are overrepresented in the dataset.
  
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
  - Display the distribution of gene expression values across different sample groups using boxplots. This feature allows users to quickly assess variability, outliers, and trends within the data.

---

## Download

1. Open your terminal.
2. Clone the repository to your local machine by running:
    ```
    git clone https://github.com/francesco-fortunato/GEVIS.git
    ```
3. Navigate to the "GEVIS" directory using the command line.
    ```
    cd GEVIS
    ```
4.  Run the following command
    ```
    docker-compose up --build
    ```
5. Go on localhost in your web browser to access the GEVIS dashboard. Utilize the interactive sliders and graphs to explore gene expression data and identify differentially expressed genes.

---

## Usage

1. **Data Upload:**
   - Upload raw gene expression data along with associated metadata, or select from available example datasets for analysis.
   
2. **Select Analysis Type:**
   - Choose between t-test, Limma, or DESeq for differential expression analysis.

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

## Docs

- [Project Report](docs/pdf/GEVIS.pdf)  
  This document contains the pdf report for our project, including methodologies, results, and conclusions.

- [Project Presentation](https://docs.google.com/presentation/d/1ZU2Z-0I1FXjg4zcnIxp9G5TC15fB4vIK7RZeQ4DTJkc/edit?usp=sharing)  
  This document contains the slides used for the project presentation.

## Contributors

- [Francesco Fortunato](https://github.com/francesco-fortunato)
- [Cristian Santaroni](https://github.com/Cristian-Santaroni)

## License

This project is licensed under the [MIT License](LICENSE).
