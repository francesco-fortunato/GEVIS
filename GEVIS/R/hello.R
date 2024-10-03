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

  data1 <- as.data.frame(data)

  # Identify the position of the 'gene' column
  gene_col_position <- which(colnames(data1) == "gene")

  # Check if 'gene' is the first column
  if (gene_col_position == 1) {
    # 'gene' is already in the first position, just set row names
    rownames(data1) <- data1[, 1]
    data1 <- data1[, -1]
  } else {
    # Move the 'gene' column to the first position
    data1 <- data1[, c(gene_col_position, setdiff(seq_along(data1), gene_col_position))]
    # Set the first column as row names
    rownames(data1) <- data1[, 1]
    # Remove the first column (it's now the row names)
    data1 <- data1[, -1]
  }

  # Convert all columns to numeric
  data1[] <- lapply(data1, function(x) as.numeric(as.character(x)))

  # Check for any non-numeric values
  if (!all(sapply(data1, is.numeric))) {
    stop("Data contains non-numeric values.")
  }

  # Print column names and data types of data1 for debugging
  print("Column names of data1:")
  print(colnames(data1))
  print("Data types of data1 columns:")
  print(sapply(data1, class))

  gsmC <- colnames(dataC)
  gsmN <- colnames(dataN)

  # Convert GSM headers to vectors and remove 'gene'
  gsmC_vector <- setdiff(gsmC, "gene")
  gsmN_vector <- setdiff(gsmN, "gene")

  # Print GSM vectors for debugging
  print("GSM C Vector:")
  print(gsmC_vector)
  print("GSM N Vector:")
  print(gsmN_vector)

  # Check if the columns in gsmC_vector and gsmN_vector exist in data1
  missing_cols_C <- setdiff(gsmC_vector, colnames(data1))
  missing_cols_N <- setdiff(gsmN_vector, colnames(data1))

  if (length(missing_cols_C) > 0 || length(missing_cols_N) > 0) {
    stop("Some columns from GSM vectors are missing in data1.")
  }

  # Combine data from case and normal samples
  data1 <- t(data1[, c(gsmC_vector, gsmN_vector)])

  # Check for any non-numeric values after transposing
  if (!all(sapply(data1, is.numeric))) {
    stop("Data contains non-numeric values after transposing.")
  }

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
  scores_var = as.data.frame(scores_var)

  # Get gene names
  gene_names <- rownames(scores_var)
  # print(gene_names)

  # Calculate the proportion of variance explained by each principal component
  explained_variance <- summary(pca)$importance[2, ] * 100  # Multiply by 100 to get percentage


  # Return the explained variance along with scores_df and scores_var
  return(list(scores_df = scores_df, scores_var = scores_var, explained_variance = explained_variance))
}



enrichment <- function(data, direction) {
  library(forcats)
  library(stringr)
  library(enrichR)

  top_term <- 10
  thr_pval <- 0.05
  dbs <- c("DisGeNET", "GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")

  list <- split(data$gene, data$direction)

  df <- lapply(list, function(x) {
    enrichr(x, dbs)
  })

  DisGeNET <- df[[direction]]$DisGeNET
  DisGeNET <- DisGeNET[DisGeNET$Adjusted.P.value < thr_pval, ]
  DisGeNET <- DisGeNET[order(DisGeNET$Adjusted.P.value, decreasing = FALSE),]

  # Convert Adjusted.P.value to string
  DisGeNET$Adjusted.P.value <- as.character(DisGeNET$Adjusted.P.value)

  DisGeNET$Gene_count <- sapply(DisGeNET$Genes, function(x) {
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  DisGeNET$Gene_ratio <- unlist(lapply(DisGeNET$Overlap, function(x) {
    total <- as.numeric(strsplit(x, "/")[[1]][2])
    count <- as.numeric(strsplit(x, "/")[[1]][1])
    Gene_ratio <- count / total
  }))

  if (length(top_term) != 0 & top_term <= nrow(DisGeNET)) {
    annotation_top <- DisGeNET[1:top_term,]
  } else {
    annotation_top <- DisGeNET
  }

  # Process GO Biological Process
  BP <- df[[direction]]$GO_Biological_Process_2021
  BP$Term <- gsub("\\s*\\(GO:\\d+\\)$", "", BP$Term)
  BP <- BP[BP$Adjusted.P.value < thr_pval, ]
  BP <- BP[order(BP$Adjusted.P.value, decreasing = FALSE),]

  # Convert Adjusted.P.value to string
  BP$Adjusted.P.value <- as.character(BP$Adjusted.P.value)

  BP$Gene_count <- sapply(BP$Genes, function(x) {
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  BP$Gene_ratio <- unlist(lapply(BP$Overlap, function(x) {
    total <- as.numeric(strsplit(x, "/")[[1]][2])
    count <- as.numeric(strsplit(x, "/")[[1]][1])
    Gene_ratio <- count / total
  }))

  if (length(top_term) != 0 & top_term <= nrow(BP)) {
    annotation_top1 <- BP[1:top_term,]
  } else {
    annotation_top1 <- BP
  }

  # Process GO Molecular Function
  MF <- df[[direction]]$GO_Molecular_Function_2021
  MF$Term <- gsub("\\s*\\(GO:\\d+\\)$", "", MF$Term)
  MF <- MF[MF$Adjusted.P.value < thr_pval, ]
  MF <- MF[order(MF$Adjusted.P.value, decreasing = FALSE),]

  # Convert Adjusted.P.value to string
  MF$Adjusted.P.value <- as.character(MF$Adjusted.P.value)

  MF$Gene_count <- sapply(MF$Genes, function(x) {
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  MF$Gene_ratio <- unlist(lapply(MF$Overlap, function(x) {
    total <- as.numeric(strsplit(x, "/")[[1]][2])
    count <- as.numeric(strsplit(x, "/")[[1]][1])
    Gene_ratio <- count / total
  }))

  if (length(top_term) != 0 & top_term <= nrow(MF)) {
    annotation_top2 <- MF[1:top_term,]
  } else {
    annotation_top2 <- MF
  }

  # Process KEGG
  KEGG <- df[[direction]]$KEGG_2021_Human
  KEGG$Term <- gsub("\\s*\\(KEGG:\\d+\\)$", "", KEGG$Term)
  KEGG <- KEGG[KEGG$Adjusted.P.value < thr_pval, ]
  KEGG <- KEGG[order(KEGG$Adjusted.P.value, decreasing = FALSE),]

  # Convert Adjusted.P.value to string
  KEGG$Adjusted.P.value <- as.character(KEGG$Adjusted.P.value)

  KEGG$Gene_count <- sapply(KEGG$Genes, function(x) {
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
  })

  KEGG$Gene_ratio <- unlist(lapply(KEGG$Overlap, function(x) {
    total <- as.numeric(strsplit(x, "/")[[1]][2])
    count <- as.numeric(strsplit(x, "/")[[1]][1])
    Gene_ratio <- count / total
  }))

  if (length(top_term) != 0 & top_term <= nrow(KEGG)) {
    annotation_top3 <- KEGG[1:top_term,]
  } else {
    annotation_top3 <- KEGG
  }

  return(list(annotation_top = annotation_top, annotation_top1 = annotation_top1, annotation_top2 = annotation_top2, annotation_top3 = annotation_top3))
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
  # Load DESeq2 library
  library(DESeq2)

  # Convert the data from list to data frames if they aren't already
  dataC <- as.data.frame(dataC)
  dataN <- as.data.frame(dataN)

  # Remove the 'gene' column before combining and ensure numeric values
  dataC_numeric <- as.matrix(dataC[,-1])  # Convert to matrix and exclude 'gene' column
  dataN_numeric <- as.matrix(dataN[,-1])  # Convert to matrix and exclude 'gene' column

  # Ensure all values are numeric
  dataC_numeric <- apply(dataC_numeric, 2, as.numeric)
  dataN_numeric <- apply(dataN_numeric, 2, as.numeric)

  # Combine the numeric matrices
  data_combined <- cbind(dataC_numeric, dataN_numeric)

  # Set row names from the 'gene' column in the original data frames
  rownames(data_combined) <- dataC$gene  # Assuming 'gene' column is the same in both dataC and dataN

  # Create metadata
  coldata <- data.frame(
    condition = factor(c(rep("Case", ncol(dataC_numeric)), rep("Normal", ncol(dataN_numeric))))
  )
  rownames(coldata) <- colnames(data_combined)

  # Create DESeq2 dataset object
  dds <- DESeqDataSetFromMatrix(countData = data_combined, colData = coldata, design = ~ condition)

  # Run DESeq2 analysis
  dds <- DESeq(dds)

  # Extract results, including log2 fold change, p-values, and adjusted p-values
  res <- results(dds, contrast = c("condition", "Case", "Normal"))

  # Convert the results to a data frame
  res_df <- as.data.frame(res)

  # Add gene names to the results
  res_df$Gene <- rownames(res_df)

  # Select relevant columns for output (log2 fold change, standard error, p-value, adjusted p-value)
  res_df <- res_df[, c("Gene", "log2FoldChange", "lfcSE", "pvalue", "padj")]

  # Rename columns for clarity
  colnames(res_df) <- c("Gene", "logFC", "lfcSE", "pvalue", "pval_adj")
  res_df$pval_adj <- as.character(res_df$pval_adj)
  res_df$logFC <- as.character(res_df$logFC)

  return(res_df)
}

surv <- function (metadata,dataC, gene){
  library(survminer)   # Load survminer first
  library(survival)    # Then load survival

  print(metadata)
  print(dataC)

  # Remove the last element from gsm_cods
  gsm_cods <- names(dataC)

  # Remove the last element from gene_expression
  gene_expression <- unlist(dataC, use.names = FALSE)

  df <- data.frame(case_id = gsm_cods, counts = gene_expression)

  print(df)

  if (nrow(df) > 1) {
    df <- df[-nrow(df), ]
  }
  df$counts <- as.numeric(df$counts)

  print(df)

  medianValue = median(df$counts)

  print(medianValue)

  df$strata = ifelse (df$counts >= medianValue,"HIGH","LOW")

  df$event = ifelse(df$case_id == metadata$GSM,  metadata$Event[match(df$case_id, metadata$GSM)], "NA")

  df$time = ifelse(df$case_id == metadata$GSM, metadata$Time_to_followUp[match(df$case_id, metadata$GSM)], "NA")


  print(df)

  df$counts <-(as.numeric(df$counts))
  metadata$Event <- as.numeric(as.character(metadata$Event))
  metadata$Time_to_followUp <- as.numeric(as.character(metadata$Time_to_followUp))

  newDF = data.frame(time = as.numeric(df$time), event = as.numeric(df$event), gene = df$strata)

  print(newDF)

  fit=survfit(Surv(time,event) ~ gene, data = newDF)

  members <- c("time", "n.risk", "n.event","n.censor","surv","strata")

  last = list(unclass(fit)[members])

  ####### UTILE PER CALCOLARE IL PVALUE DEL GENE PASSATO#######
  obj_pval = surv_pvalue(fit,newDF)
  members2 = c("pval")
  pval=list(unclass(obj_pval)[members2])
  print(pval)
  #############################################################

  return(list (obj = last, pval = pval))

}
