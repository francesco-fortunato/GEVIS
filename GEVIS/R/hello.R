# Function to get current timestamp in seconds
current_time <- function() {
  as.numeric(Sys.time())
}

# Function to hash rawdata for unique key
hash <- function(data) {
  digest::digest(data)
}

# Define a cache key function
cacheKey <- function(func_name, data) {
  data_hash <- hash(data)
  key <- list(func = func_name, data_hash = data_hash)

  return(key)
}

# Function to check if cache is expired
is_cache_expired <- function(timestamp, max_age) {
  (current_time() - timestamp) > max_age
}

pval <- function(data, N, M) {
  library(R.cache)
  library(digest)

  # Generate a cache key
  key <- cacheKey("pval", list(data, N, M))

  # Check if the result is already cached
  cached_result <- loadCache(key)

  # If cache exists, return it
  if (!is.null(cached_result)) {
    print("CACHE FOUND")
    # If cache is still valid, return cached data
    return(cached_result$data)
  }

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

  # Save cache with a timestamp
  saveCache(list(timestamp = current_time(), data = result), key)

  return(result)
}


pca <- function(data, dataC, dataN) {
  library(FactoMineR)
  library(factoextra)
  library(dplyr)
  library(R.cache)
  library(digest)

  # Generate a cache key
  key <- cacheKey("pca", list(data, dataC, dataN))

  # Check if the result is already cached
  cached_result <- loadCache(key)

  # If cache exists, return it
  if (!is.null(cached_result)) {
    print("CACHE FOUND")
    # If cache is still valid, return cached data
    return(cached_result$data)
  }

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

  # Save cache with a timestamp
  saveCache(list(timestamp = current_time(), data = list(scores_df = scores_df, scores_var = scores_var, explained_variance = explained_variance)), key)


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
  library(R.cache)
  library(digest)

  # Generate a cache key
  key <- cacheKey("variation", rawdata)

  # Check if the result is already cached
  cached_result <- loadCache(key)

  # If cache exists, return it
  if (!is.null(cached_result)) {
    print("CACHE FOUND")
    # If cache is still valid, return cached data
    return(cached_result$data)
  }

  # Extract gene names
  genes <- rawdata[, ncol(rawdata)]

  # Calculate variation using IQR for each gene
  variation <- apply(rawdata[, -ncol(rawdata)], 1, IQR)

  # Create a data frame containing gene names and variation
  variat_data <- data.frame(Gene = genes, Variation = variation)

  # Convert data frame to JSON
  json_data <- toJSON(variat_data)

  # Save cache with a timestamp
  saveCache(list(timestamp = current_time(), data = json_data), key)

  # Return JSON data
  return(json_data)
}

limmaDE <- function(dataC, dataN) {
  library(limma)

  library(R.cache)
  library(digest)

  # Generate a cache key
  key <- cacheKey("deseq2DE", list(dataC, dataN))

  # Check if the result is already cached
  cached_result <- loadCache(key)

  # If cache exists, return it
  if (!is.null(cached_result)) {
    print("CACHE FOUND")
    # If cache is still valid, return cached data
    return(cached_result$data)
  }


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

  # Save cache with a timestamp
  saveCache(list(timestamp = current_time(), data = results), key)

  return(results)
}




deseq2DE <- function(dataC, dataN) {
  # Load DESeq2 library
  library(DESeq2)

  library(R.cache)
  library(digest)

  # Generate a cache key
  key <- cacheKey("deseq2DE", list(dataC, dataN))

  # Check if the result is already cached
  cached_result <- loadCache(key)

  # If cache exists, return it
  if (!is.null(cached_result)) {
    print("CACHE FOUND")
    # If cache is still valid, return cached data
    return(cached_result$data)
  }

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

  # Save cache with a timestamp
  saveCache(list(timestamp = current_time(), data = res_df), key)

  return(res_df)
}

surv <- function (metadata, dataC, gene) {
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

  # Calculate the median and convert it to character
  medianValue <- median(df$counts)
  medianValue_char <- as.character(medianValue)

  print(medianValue)

  df$strata = ifelse(df$counts >= medianValue, "HIGH", "LOW")

  df$event = ifelse(df$case_id == metadata$GSM, metadata$Event[match(df$case_id, metadata$GSM)], "NA")

  df$time = ifelse(df$case_id == metadata$GSM, metadata$Time_to_followUp[match(df$case_id, metadata$GSM)], "NA")

  print(df)

  df$counts <- as.numeric(df$counts)
  metadata$Event <- as.numeric(as.character(metadata$Event))
  metadata$Time_to_followUp <- as.numeric(as.character(metadata$Time_to_followUp))

  newDF = data.frame(time = as.numeric(df$time), event = as.numeric(df$event), gene = df$strata)

  print(newDF)

  fit <- survfit(Surv(time, event) ~ gene, data = newDF)

  members <- c("time", "n.risk", "n.event", "n.censor", "surv", "strata")
  last <- list(unclass(fit)[members])

  ####### Calculate p-value for the passed gene #######
  obj_pval <- surv_pvalue(fit, newDF)
  members2 <- c("pval")
  pval <- list(unclass(obj_pval)[members2])
  print(pval)
  ####################################################

  # Return the survival object, p-value, and median value (as character)
  return(list(obj = last, pval = pval, median = medianValue_char))
}


heatmap <- function(data, metadata, field, case, distance, method, show_cols, show_rows, cutreeCols, cutreeRows) {
  # Ensure necessary libraries are loaded
  library(pheatmap)
  library(grid)  # Important for custom text graphics
  library(gridExtra)  # For grid.arrange
  library(ggplot2)  # For ggsave
  library(RColorBrewer)  # For color palettes
  library(viridis)  # For better color scales

  print(show_cols)

  # Convert data to a data frame and ensure numeric values
  data <- as.data.frame(data)
  data[,-ncol(data)] <- lapply(data[,-ncol(data)], as.numeric)

  # Set row names from the last column and remove it
  rownames(data) <- data[, ncol(data)]
  data <- data[, -ncol(data)]

  print("Processed Data:")
  print(data)

  # Transpose the metadata for easier access
  metadata <- as.data.frame(t(metadata))

  # Set the first row as column names
  colnames(metadata) <- metadata[1, ]
  metadata <- metadata[-1, , drop = FALSE]  # Remove the first row now that it's column names

  # Check if the primary field is a valid column in metadata
  if (!field %in% colnames(metadata)) {
    stop(paste("Error: The field does not exist in metadata:", field))
  }

  # Filter metadata fields to exclude those with more than 10 unique values
  metadata_fields <- metadata[, colnames(metadata), drop = FALSE]  # Extract as a data frame
  unique_counts <- sapply(metadata_fields, function(col) length(unique(col)))  # Count unique values

  # Keep columns with <= 10 unique values and check numeric fields
  metadata_filtered <- metadata_fields[, unique_counts <= 10, drop = FALSE]

  # Check fields with > 10 unique values for numeric types and keep them for scaling
  numeric_metadata <- metadata_fields[, unique_counts > 10, drop = FALSE]

  # Identify numeric columns
  numeric_cols <- sapply(numeric_metadata, function(col) {
    is_numeric <- suppressWarnings(!any(is.na(as.numeric(as.character(col)))))  # Check if conversion to numeric is possible
    return(is_numeric)
  })

  # Ensure numeric_cols is a logical vector
  if (length(numeric_cols) == 0) {
    numeric_cols <- logical(0)  # Handle case where there are no columns
  } else {
    numeric_cols <- as.logical(numeric_cols)  # Convert to logical vector
  }

  # Only keep numeric columns with > 10 unique values
  if (any(numeric_cols)) {
    numeric_metadata_filtered <- numeric_metadata[, numeric_cols, drop = FALSE]
  } else {
    numeric_metadata_filtered <- data.frame()  # Create an empty data frame if no numeric columns
  }

  # Check if we have any numeric columns to include
  if (ncol(numeric_metadata_filtered) > 0) {
    metadata_final <- cbind(metadata_filtered, numeric_metadata_filtered)
  } else {
    metadata_final <- metadata_filtered
  }

  # Create annotation colors dynamically for categorical fields
  annotation_colors <- list()
  for (name in names(metadata_final)) {
    if (name %in% names(numeric_metadata_filtered)) {
      next
    }

    if (name == field) {
      levels <- unique(metadata_final[[name]])
      if (case == levels[1]) {
        annotation_colors[[name]] <- c(RColorBrewer::brewer.pal(3, "Set1")[1], RColorBrewer::brewer.pal(3, "Set1")[3])
      } else {
        annotation_colors[[name]] <- c(RColorBrewer::brewer.pal(3, "Set1")[3], RColorBrewer::brewer.pal(3, "Set1")[1])
      }
      names(annotation_colors[[name]]) <- levels
    } else {
      levels <- unique(metadata_final[[name]])
      annotation_colors[[name]] <- RColorBrewer::brewer.pal(length(levels), "Set3")[1:length(levels)]
      names(annotation_colors[[name]]) <- levels
    }
  }

  if (ncol(numeric_metadata_filtered) > 0) {
    for (name in names(numeric_metadata_filtered)) {
      annotation_colors[[name]] <- colorRampPalette(c("lightblue", "blue"))(100)
    }
    metadata_final[names(numeric_metadata_filtered)] <- lapply(metadata_final[names(numeric_metadata_filtered)], as.numeric)
  }

  # Function to truncate row names
  truncate_rownames <- function(rownames, max_length) {
    sapply(rownames, function(name) {
      if (nchar(name) > max_length) {
        return(paste0(substr(name, 1, max_length), "..."))
      } else {
        return(name)
      }
    })
  }

  max_length <- 15
  truncated_rownames <- truncate_rownames(rownames(data), max_length)
  truncated_rownames <- make.unique(truncated_rownames)
  rownames(data) <- truncated_rownames

  # Determine font size based on number of columns
  n_cols <- ncol(data)
  col_font_size <- ifelse(n_cols <= 25, 12,   # If 5 or fewer columns, use font size 12
                          ifelse(n_cols <= 50, 8,  # If between 26 and 50, use font size 10
                                 ifelse(n_cols <= 75, 6,  # If between 26 and 50, use font size 10
                                        ifelse(n_cols <= 120, 4,  # If between 26 and 100, use font size 10
                                    2))))  # For more than 10 columns, use font size 8

  # Custom function to draw diagonal column names with smaller font size
  draw_colnames_45 <- function(coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size

    # Set character size (cex) to make the text smaller
    res = textGrob(coln, x = x,
                   y = unit(1, "npc") - unit(4, "bigpts"),  # Adjust y position if needed
                   vjust = 0.5, hjust = 1,
                   rot = 90,
                   gp = gpar(..., lineheight = 1, cex = 1, fontsize = col_font_size, fontfamily = "Arial"))  # Use dynamic font size
    return(res)
  }

  n_rows <- nrow(data)
  row_font_size <- ifelse(n_rows <= 25, 12,   # If 5 or fewer columns, use font size 12
                          ifelse(n_rows <= 50, 8,  # If between 26 and 50, use font size 10
                                 ifelse(n_rows <= 120, 5,  # If between 26 and 100, use font size 10
                                        ifelse(n_rows <= 170, 4,  # If between 26 and 100, use font size 10
                                               ifelse(n_rows <= 250, 2.6,  # If between 26 and 100, use font size 10
                                                      ifelse(n_rows <= 350, 2,  # If between 26 and 100, use font size 10
                                                             1.5))))))  # For more than 10 columns, use font size 8


  # Custom function to draw row names with dynamic font size
  draw_rownames_custom <- function(rown, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(rown), gaps)
    y = rev(coord$coord) - 0.5 * rev(coord$size)  # Reverse for correct alignment

    # Set character size (cex) to make the text smaller
    res = textGrob(rown, x = unit(1, "npc") - unit(4, "bigpts"),  # Adjust x position if needed
                   y = y,
                   vjust = 0.5, hjust = 1,
                   gp = gpar(..., lineheight = 1, cex = 1, fontsize = row_font_size, fontfamily = "Arial"))  # Use dynamic font size
    return(res)
  }


  # Overwrite default draw_colnames with the custom version
  assignInNamespace(x = "draw_colnames", value = draw_colnames_45, ns = asNamespace("pheatmap"))

  # Generate the pheatmap plot
  my_pheatmap <- pheatmap(
    data,
    scale = "row",
    border_color = NA,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = distance,
    clustering_distance_cols = distance,
    clustering_method = method,
    annotation_col = metadata_final,
    color = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(1000),
    show_rownames = show_rows,
    show_colnames = show_cols,
    cutree_cols = cutreeCols,
    cutree_rows = cutreeRows,
    annotation_colors = annotation_colors,
    fontsize_row = row_font_size
  )

  # Save the plot as SVG using ggsave
  ggsave(filename = "heatmap_output.svg", plot = my_pheatmap[[4]], width = 10, height = 10, dpi=72)
  rm(data)
  rm(metadata)
  return(TRUE)
}


heatmap_mod <- function(data, metadata, field, case, distance, method, show_cols, show_rows, cutreeCols, cutreeRows) {
  # Load necessary libraries
  library(pheatmap)
  library(grid)  # Important for custom text graphics
  library(gridExtra)  # For grid.arrange
  library(ggplot2)  # For ggsave
  library(RColorBrewer)  # For color palettes
  library(viridis)  # For better color scales
  library(jsonlite)

  # Generate pheatmap to obtain the clustering information
  #pheat <- pheatmap(data_matrix)
  # Convert data to a data frame and ensure numeric values
  data <- as.data.frame(data)
  data[,-ncol(data)] <- lapply(data[,-ncol(data)], as.numeric)

  # Set row names from the last column and remove it
  rownames(data) <- data[, ncol(data)]
  data <- data[, -ncol(data)]

  print("Processed Data:")
  print(data)

  # Transpose the metadata for easier access
  metadata <- as.data.frame(t(metadata))

  # Set the first row as column names
  colnames(metadata) <- metadata[1, ]
  metadata <- metadata[-1, , drop = FALSE]  # Remove the first row now that it's column names

  # Check if the primary field is a valid column in metadata
  if (!field %in% colnames(metadata)) {
    stop(paste("Error: The field does not exist in metadata:", field))
  }

  # Filter metadata fields to exclude those with more than 10 unique values
  metadata_fields <- metadata[, colnames(metadata), drop = FALSE]  # Extract as a data frame
  unique_counts <- sapply(metadata_fields, function(col) length(unique(col)))  # Count unique values

  # Keep columns with <= 10 unique values and check numeric fields
  metadata_filtered <- metadata_fields[, unique_counts <= 10, drop = FALSE]

  # Check fields with > 10 unique values for numeric types and keep them for scaling
  numeric_metadata <- metadata_fields[, unique_counts > 10, drop = FALSE]

  # Identify numeric columns
  numeric_cols <- sapply(numeric_metadata, function(col) {
    is_numeric <- suppressWarnings(!any(is.na(as.numeric(as.character(col)))))  # Check if conversion to numeric is possible
    return(is_numeric)
  })

  # Ensure numeric_cols is a logical vector
  if (length(numeric_cols) == 0) {
    numeric_cols <- logical(0)  # Handle case where there are no columns
  } else {
    numeric_cols <- as.logical(numeric_cols)  # Convert to logical vector
  }

  # Only keep numeric columns with > 10 unique values
  if (any(numeric_cols)) {
    numeric_metadata_filtered <- numeric_metadata[, numeric_cols, drop = FALSE]
  } else {
    numeric_metadata_filtered <- data.frame()  # Create an empty data frame if no numeric columns
  }

  # Check if we have any numeric columns to include
  if (ncol(numeric_metadata_filtered) > 0) {
    metadata_final <- cbind(metadata_filtered, numeric_metadata_filtered)
  } else {
    metadata_final <- metadata_filtered
  }

  # Create annotation colors dynamically for categorical fields
  annotation_colors <- list()
  for (name in names(metadata_final)) {
    if (name %in% names(numeric_metadata_filtered)) {
      next
    }

    if (name == field) {
      levels <- unique(metadata_final[[name]])
      if (case == levels[1]) {
        annotation_colors[[name]] <- c(RColorBrewer::brewer.pal(3, "Set1")[1], RColorBrewer::brewer.pal(3, "Set1")[3])
      } else {
        annotation_colors[[name]] <- c(RColorBrewer::brewer.pal(3, "Set1")[3], RColorBrewer::brewer.pal(3, "Set1")[1])
      }
      names(annotation_colors[[name]]) <- levels
    } else {
      levels <- unique(metadata_final[[name]])
      annotation_colors[[name]] <- RColorBrewer::brewer.pal(length(levels), "Set3")[1:length(levels)]
      names(annotation_colors[[name]]) <- levels
    }
  }

  if (ncol(numeric_metadata_filtered) > 0) {
    for (name in names(numeric_metadata_filtered)) {
      annotation_colors[[name]] <- colorRampPalette(c("lightblue", "blue"))(100)
    }
    metadata_final[names(numeric_metadata_filtered)] <- lapply(metadata_final[names(numeric_metadata_filtered)], as.numeric)
  }

  pheat <- pheatmap(
    data,
    scale = "row",
    border_color = NA,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = distance,
    clustering_distance_cols = distance,
    clustering_method = method,
    annotation_col = metadata_final,
    color = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(1000),
    show_rownames = show_rows,
    show_colnames = show_cols,
    cutree_cols = cutreeCols,
    cutree_rows = cutreeRows,
    annotation_colors = annotation_colors
  )

  # Extract row and column clustering
  row_clust <- pheat$tree_row
  col_clust <- pheat$tree_col

  # Helper function to convert dendrogram to a nested list
  dendrogram_to_list <- function(d) {
    if (is.leaf(d)) {
      return(list(name = attr(d, "label")))
    } else {
      return(list(
        children = lapply(d, dendrogram_to_list),
        height = attr(d, "height")
      ))
    }
  }

  # Convert row and column dendrograms to lists
  row_dendro_list <- dendrogram_to_list(as.dendrogram(row_clust))
  col_dendro_list <- dendrogram_to_list(as.dendrogram(col_clust))

  # Prepare the cell data in a JSON-compatible format
  cell_data_list <- list()
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      cell_data_list[[length(cell_data_list) + 1]] <- list(
        row = i,
        col = j,
        value = data[i, j]
      )
    }
  }

  # Extract cluster assignments
  row_clusters <- cutree(row_clust, k = cutreeRows)
  col_clusters <- cutree(col_clust, k = cutreeCols)

  # Convert to data frames
  row_clusters_df <- data.frame(Row = names(row_clusters), Cluster = row_clusters)
  col_clusters_df <- data.frame(Col = names(col_clusters), Cluster = col_clusters)

  # Combine all information into a single list
  heatmap_json <- list(
    row_dendrogram = row_dendro_list,
    col_dendrogram = col_dendro_list,
    cell_data = cell_data_list,
    row_clusters = row_clusters_df,
    col_clusters = col_clusters_df
  )

  # Convert the combined list to JSON
  heatmap_json_str <- toJSON(heatmap_json)

  # Save the JSON string to a file
  return(heatmap_json_str)
}
