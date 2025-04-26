calculate_difference <- function(dataC, dataN) {
  rowmeans_diff <- rowMeans(dataC) - rowMeans(dataN)
  return(rowmeans_diff)
}
