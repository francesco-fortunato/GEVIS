hello <- function(data) {
  # Save the column of genes before removing the last column
  genes <- data[, ncol(data)]

  # Assuming data, dataN, and dataC are your matrices
  data <- data[, -ncol(data)]


  N <- 49
  M <- 58
  # Now, the last column has been removed from each matrix

  all_equal <- function(x) {
    return(length(unique(x)) == 1)
  }

  # Compute p-values with checks for constant data
  pval <- apply(data, 1, function(x) {
    if (all_equal(x[1:N]) || all_equal(x[(N+1):(M+N)])) {
      return(1)
    } else {
      return(t.test(x[1:N], x[(N+1):(M+N)], paired = FALSE)$p.value)
    }
  })

  # Adjustment p-value
  pval_adj <- p.adjust(pval, method="fdr")

  # Add the column of genes back to the result
  result <- data.frame(Gene = genes, pval_adj = pval_adj)

  return(result)
}


pval <- function(data, N, M) {
  # Save the column of genes before removing the last column
  genes <- data[, ncol(data)]

  # Assuming data, dataN, and dataC are your matrices
  data <- data[, -ncol(data)]

  # Define a helper function to check if all elements are equal
  all_equal <- function(x) {
    return(length(unique(x)) == 1)
  }

  # Define a helper function to check for low variance
  low_variance <- function(x, threshold = 1e-10) {
    return(var(x) < threshold)
  }

  # Compute p-values with checks for constant data and low variance
  pval <- apply(data, 1, function(x) {
    if (all_equal(x[1:N]) || all_equal(x[(N+1):(M+N)]) || all_equal(x[1:(M+N)])) {
      return(1)
    } else {
      return(t.test(x[1:N], x[(N+1):(M+N)], paired = FALSE)$p.value)
    }
  })

  # Adjustment p-value
  pval_adj <- p.adjust(pval, method="fdr")

  # Convert pval_adj to character (string)
  pval_adj_str <- as.character(pval_adj)

  # Add the column of genes back to the result
  result <- data.frame(Gene = genes, pval_adj = pval_adj_str)

  return(result)
}


pca <- function(data, dataC, dataN) {
  library(FactoMineR)
  library(factoextra)
  library(dplyr)

  data1 <- as.data.frame.list(data)

  # Move the last column to the first position
  data1 <- data1[, c(ncol(data1), 1:(ncol(data1)-1))]

  # Set the first column as row names
  rownames(data1) <- data1[, 1]

  # Remove the first column (it's now the row names)
  data1 <- data1[, -1]

  gsmC <- colnames(dataC)
  gsmN <- colnames(dataN)

  # Convert GSM headers to vectors and remove the last element
  gsmC_vector <- head(unlist(gsmC), -1)
  gsmN_vector <- head(unlist(gsmN), -1)

  # Combine data from case and normal samples
  data1 <- t(data1[, c(gsmC_vector, gsmN_vector)])
  # Create groups vector
  groups <- c(rep("case", length(gsmC_vector)), rep("normal", length(gsmN_vector)))

  # Perform PCA
  pca <- prcomp(data1, center = TRUE, scale. = TRUE, retx = TRUE)

  # Compute scores
  scores <- pca$x

  # Convert scores to a data frame
  scores_df <- as.data.frame(scores)
  scores_df$Group <- groups

  # Get PCA variable contributions
  scores_var <- get_pca_var(pca)$contrib
  colnames(scores_var) <- paste0('PC', seq(1, ncol(scores_var)))

  # Sort by PC1 contribution
  scores_var <- scores_var[order(scores_var[,"PC1"], decreasing = TRUE), ]
  scores_var= as.data.frame(scores_var)

  # Get gene names
  gene_names <- rownames(scores_var)
  # print(gene_names)

  # Return a list containing scores_df and scores_var
  return(list(scores_df = scores_df, scores_var = scores_var))
}

enrichment <- function (data,direction){
  library(forcats)
  library(stringr)
  library(enrichR)

  top_term <- 10
  thr_pval <- 0.05
  dbs <- c("DisGeNET","GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")

  list <- split(data$gene,data$direction)

  df <- lapply(list, function(x){
    enrichr(x, dbs)
  })



  DisGeNET <- df[[direction]]$DisGeNET
  DisGeNET <- DisGeNET[DisGeNET$Adjusted.P.value < thr_pval, ]
  DisGeNET <- DisGeNET[order(DisGeNET$Adjusted.P.value, decreasing = F),]

  DisGeNET$Gene_count <- sapply(DisGeNET$Genes, function(x){
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  DisGeNET$Gene_ratio <- unlist(lapply(DisGeNET$Overlap, function(x){

    total <- as.numeric(strsplit(x,"/")[[1]][2])
    count <- as.numeric(strsplit(x,"/")[[1]][1])

    Gene_ratio <- count/total

  } ))

  if(length(top_term) != 0 & top_term <= nrow(DisGeNET)){
    annotation_top <- DisGeNET[1:top_term,]
  }else{
    annotation_top <- DisGeNET
  }



  BP <- df[[direction]]$GO_Biological_Process_2021
  BP$Term <- gsub("\\s*\\(GO:\\d+\\)$", "", BP$Term)
  BP <- BP[BP$Adjusted.P.value < thr_pval, ]
  BP <- BP[order(BP$Adjusted.P.value, decreasing = F),]

  BP$Gene_count <- sapply(BP$Genes, function(x){
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  BP$Gene_ratio <- unlist(lapply(BP$Overlap, function(x){

    total <- as.numeric(strsplit(x,"/")[[1]][2])
    count <- as.numeric(strsplit(x,"/")[[1]][1])

    Gene_ratio <- count/total

  } ))

  if(length(top_term) != 0 & top_term <= nrow(BP)){
    annotation_top1 <- BP[1:top_term,]
  }else{
    annotation_top1 <- BP
  }

  MF <- df[[direction]]$GO_Molecular_Function_2021
  MF$Term <- gsub("\\s*\\(GO:\\d+\\)$", "", MF$Term)
  MF <- MF[MF$Adjusted.P.value < thr_pval, ]
  MF <- MF[order(MF$Adjusted.P.value, decreasing = F),]

  MF$Gene_count <- sapply(MF$Genes, function(x){
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  MF$Gene_ratio <- unlist(lapply(MF$Overlap, function(x){

    total <- as.numeric(strsplit(x,"/")[[1]][2])
    count <- as.numeric(strsplit(x,"/")[[1]][1])

    Gene_ratio <- count/total

  } ))

  if(length(top_term) != 0 & top_term <= nrow(MF)){
    annotation_top2 <- MF[1:top_term,]
  }else{
    annotation_top2 <- MF
  }

  KEGG <- df[[direction]]$KEGG_2021_Human
  KEGG$Term <- gsub("\\s*\\(KEGG:\\d+\\)$", "", KEGG$Term)
  KEGG <- KEGG[KEGG$Adjusted.P.value < thr_pval, ]
  KEGG <- KEGG[order(KEGG$Adjusted.P.value, decreasing = F),]

  KEGG$Gene_count <- sapply(KEGG$Genes, function(x){
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  KEGG$Gene_ratio <- unlist(lapply(KEGG$Overlap, function(x){

    total <- as.numeric(strsplit(x,"/")[[1]][2])
    count <- as.numeric(strsplit(x,"/")[[1]][1])

    Gene_ratio <- count/total

  } ))

  if(length(top_term) != 0 & top_term <= nrow(KEGG)){
    annotation_top3 <- KEGG[1:top_term,]
  }else{
    annotation_top3 <- KEGG
  }

  return(list(annotation_top = annotation_top, annotation_top1 = annotation_top1, annotation_top2= annotation_top2,annotation_top3= annotation_top3))

}


variation <- function(rawdata) {
  library(jsonlite)

  # Extract gene names
  genes <- rawdata[, ncol(rawdata)]

  # Calculate variation using IQR for each gene
  variation <- apply(rawdata[, -ncol(rawdata)], 1, IQR)

  # Create a data frame containing gene names and variation
  variat_data <- data.frame(Gene = genes, Variation = variation)

  # Convert data frame to JSON
  json_data <- toJSON(variat_data)

  # Return JSON data
  return(json_data)
}

limmaDE <- function(dataC, dataN) {
  library(limma)

  # Convert lists to data frames
  dataC <- as.data.frame(dataC)
  dataN <- as.data.frame(dataN)

  # Convert all columns to numeric except for the last column (which is assumed to be the gene column)
  dataC[,-ncol(dataC)] <- lapply(dataC[,-ncol(dataC)], as.numeric)
  dataN[,-ncol(dataN)] <- lapply(dataN[,-ncol(dataN)], as.numeric)

  # Set row names to gene names (assuming the gene names are in the last column)
  rownames(dataC) <- dataC[,ncol(dataC)]
  rownames(dataN) <- dataN[,ncol(dataN)]

  # Drop the gene column
  dataC <- dataC[,-ncol(dataC)]
  dataN <- dataN[,-ncol(dataN)]

  # Combine data
  data_combined <- cbind(dataC, dataN)

  # Create group vector
  group <- factor(c(rep("Case", ncol(dataC)), rep("Normal", ncol(dataN))))

  # Design matrix
  design <- model.matrix(~ group)

  # Fit the model
  fit <- lmFit(data_combined, design)
  fit <- eBayes(fit)

  # Get the top table
  results <- topTable(fit, coef = 2, number = Inf)

  # Rename column
  colnames(results)[which(colnames(results) == "adj.P.Val")] <- "pval_adj"

  # Convert pval_adj column to character
  results$pval_adj <- as.character(results$pval_adj)

  return(results)
}




deseq2DE <- function(dataC, dataN) {
  library(DESeq2)

  # Convert lists to data frames
  dataC <- as.data.frame(dataC)
  dataN <- as.data.frame(dataN)

  # Convert all columns to numeric except for the last column (which is assumed to be the gene column)
  dataC[,-ncol(dataC)] <- lapply(dataC[,-ncol(dataC)], as.numeric)
  dataN[,-ncol(dataN)] <- lapply(dataN[,-ncol(dataN)], as.numeric)

  # Set row names to gene names (assuming the gene names are in the last column)
  rownames(dataC) <- dataC[,ncol(dataC)]
  rownames(dataN) <- dataN[,ncol(dataN)]

  # Drop the gene column
  dataC <- dataC[,-ncol(dataC)]
  dataN <- dataN[,-ncol(dataN)]

  # Combine data
  data_combined <- cbind(dataC, dataN)

  # Create metadata
  coldata <- data.frame(
    condition = factor(c(rep("Case", ncol(dataC)), rep("Normal", ncol(dataN))))
  )
  rownames(coldata) <- colnames(data_combined)

  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = data_combined, colData = coldata, design = ~ condition)

  # Run DESeq2
  dds <- DESeq(dds)

  # Get results
  res <- results(dds, contrast = c("condition", "Case", "Normal"))

  # Adjust p-values (This step is not necessary as DESeq2 automatically calculates adjusted p-values)
  # res$padj <- p.adjust(res$pvalue, method = "BH")

  # Convert to data frame
  res_df <- as.data.frame(res)

  # Add gene names
  res_df$Gene <- rownames(res_df)

  # Select relevant columns
  res_df <- res_df[, c("Gene", "log2FoldChange", "lfcSE", "pvalue", "padj")]

  # Convert pval_adj column to character
  res_df$padj <- as.character(res_df$padj)

  # Rename columns
  colnames(res_df) <- c("Gene", "logFC", "lfcSE", "pvalue", "pval_adj")

  return(res_df)
}
