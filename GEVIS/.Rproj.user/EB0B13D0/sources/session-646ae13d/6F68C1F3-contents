hello <- function(data, dataC, dataN, variation, thr_prc) {
  ind <- which(variation <= thr_prc)
  number <- nrow(data)
  print(paste("genes before filtering:", number))

  newdata <- data[-ind, ]
  newdataC <- dataC[-ind, ]
  newdataN <- dataN[-ind, ]

  # Exclude the first column (gene names) before calculating row means
  newdataC_numeric <- as.matrix(newdataC[, -1])
  newdataN_numeric <- as.matrix(newdataN[, -1])

  number <- nrow(newdata)

  # Calculate row means for numeric columns only
  logFC <- rowMeans(newdataC_numeric) - rowMeans(newdataN_numeric)

  print(paste("genes after filtering:", number))

  # Create a named vector with gene names as names and log fold change values as elements
  logFC_named <- setNames(logFC, newdataC[, 1])

  return(list(newdata = newdata, newdataC = newdataC, newdataN = newdataN, logFC = logFC_named))
}
