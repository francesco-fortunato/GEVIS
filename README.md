# GEVIS - Gene Expression VIS

## Description

This project implements a dashboard for conducting a differential gene expression analysis on the [GSE10072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072) GEO dataset. The dashboard provides interactive visualization and analysis tools to explore gene expression data and identify differentially expressed genes between two biological conditions.

<img src="pictures/GEVIS Dashboard.png" alt="Dashboard Image" style="max-width: 100%; height: auto;">

## Features

- **Interquartile Range (IQR) Graph:** Visualizes the interquartile range of gene expression data with a slider for selecting the desired range.
- **Log Transformation:** Applies a logarithmic transformation to the gene expression data to stabilize variance.
- **LogFC Distribution Graph:** Displays the distribution of log fold change (LogFC) values with a slider for selecting the LogFC threshold.
- **Statistical Testing:** Conducts a t-test to identify genes that show significant differences in expression levels between conditions. P-values and adjusted p-values using false discovery rate (FDR) are calculated.
- **Scatterplot:** Generates a scatterplot with each circle representing a gene, colored based on LogFC values. Upregulated genes are colored blue, downregulated genes are colored red, and genes that do not meet the LogFC or adjusted p-value thresholds are colored grey. Clicking on a gene triggers an animation and adds the corresponding axis in the scatterplot. Additionally, genes can be selected from a dropdown menu for further analysis.
- **Boxplot:** Displays the distribution of gene expression values across samples for selected genes. Users can interactively select genes from the scatterplot, and the corresponding boxplot will update to show the distribution of expression values for those genes. The boxplot allows users to visualize the variability in gene expression between sample groups and identify potential outliers or patterns.
- **Parallel Coordinate Plot:** Visualizes gene expression profiles across samples using parallel coordinates. Lines represent individual samples, and users can add or remove metadata axes, select encoding options, and interactively reorder axes. Clicking on a jitter in the scatterplot highlights the corresponding line in the parallel coordinate plot and viceversa.
- **Principal Component Analysis (PCA):** Performs PCA on the gene expression data to identify potential clusters or patterns among samples.
- **Heatmap:** Displays the variable contributions of genes to principal components computed through PCA analysis, facilitating a detailed examination of gene contributions to sample clustering and variability.

## Usage

1. Open your terminal.
2. Clone the repository to your local machine by running:
    ```
    git clone https://github.com/francesco-fortunato/GEVIS.git
    ```
3. Navigate to the "GEVIS" directory using the command line.
    ```
    cd GEVIS
    ```
4. Install dependencies by running:
    ```
    npm install
    ```
5. Open a new terminal, pull the opencpu/rstudio image and launch a Docker container named "mybox" with the OpenCPU/RStudio image using the command:
    ```
    docker pull opencpu/rstudio
    docker run --name mybox -t -p 80:80 opencpu/rstudio
    ```
6. Open another terminal window, then execute the following commands:
    ```
    docker exec -i -t mybox /bin/bash
    sudo -i
    apt-get update
    apt-get install cmake
    ```
7. Access RStudio by visiting localhost/rstudio/ (login credentials: username - opencpu, password - opencpu).
8. Create a package, calling it "GEVIS", and paste inside the hello.R function the content of the file that you will find [here](GEVIS/R/hello.R). To build the project, press Ctrl+Shift+B.
9. Back in the terminal from step 4, run the command:
    ```
    node ./server.js
    ```
10. Open localhost:3000 in your web browser to access the GEVIS dashboard. Utilize the interactive sliders and graphs to explore gene expression data and identify differentially expressed genes.

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
