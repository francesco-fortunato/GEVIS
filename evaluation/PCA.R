rm(list = ls())

library(ggbiplot)
library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)
library(survminer)

setwd("~/")
path_in <- "~/Results/"
################################################

dirRes = "Results/"
if(!dir.exists(dirRes)){
  dir.create(dirRes)
}else {
  print(paste("the directory ", dirRes,"already exist"))
}

dataset = "RA" # a folder for my case of study
dirDataset = paste0(dirRes,dataset,"/")
if(!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else {
  print(paste("the directory ", dirDataset,"already exist"))
}

#file score plot is the plot of the score of pc1 against pc2 that are the old variable in the new
# reference syste,
file_score_plot <- paste0(dirDataset, "score_plot.pdf")

#pareto scree plot is the plot about the expplained variance of each component
# how many component could be extends to retain?
file_pareto_scree_plot <- paste0(dirDataset, "pareto_scree_plot.pdf")
#loading plot is a mess

#contribution of each variable in each component
#for each component which are the gene that are contributing much? we select them
file_contribution_plot <- paste0(dirDataset, "PC_contribution_plot.pdf")

################################################
# 1. Importing data

data <- read.table(paste0(path_in, dataset, "/matrix_DEG.txt"),
                   header = T, sep = "\t", quote = "",
                   check.names = F, row.names = 1)

list_normal <- read.table (paste0(path_in, dataset, "/normal.txt"),
                           header = F, quote = "", check.names = F)$V1

list_case <- read.table (paste0(path_in, dataset, "/case.txt"),
                         header = F, quote = "", check.names = F)$V1

data <- t(data[,c(list_case, list_normal)])

groups <- c(rep("case", length(list_case)), rep("normal", length(list_normal)))
################################################
# 2. Apply PCA
# Rows of data correspond to observations (samples), columns to variables (geni)

pca <- prcomp(data, center = T, scale. = T, retx = T)
################################################
# 3. Compute score and score plot
# (scores = the coordinates of old data (observations) in the new systems, that are the PCs)

# pca$x = t(data)*pca$rotation
scores <- pca$x

# alternative for score computation
# scores <- get_pca_ind(pca)$coord
# colnames(scores) <- paste0('PC', seq(1,ncol(scores)))

pdf(file_score_plot,width=5, height=5)
g <- ggbiplot(pca,
              obs.scale = 1,
              var.axes = F,
              ellipse = T,
              groups = groups)
print(g)
dev.off()

# alternative score plot
fviz_pca_ind(pca,
             col.ind = groups, # "cos2"
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE,
             repel = TRUE, # Avoid text overlapping
)

################################################
# 4. Compute eigenvalue
# eigenvalues of the covariance matrix ordered in decreasing order (from the largest to the smallest)
eigenvalue = pca$sdev^2

# variance explained by each PC
varS <- round(eigenvalue/sum(eigenvalue)*100, 2)
names(varS) = paste0('PC', seq(1,length(varS)))

pdf(file_pareto_scree_plot, width =5, height=5)

# pareto chart
pareto.chart(varS[1:10])

# scree plot
fviz_eig(pca,addlabels = TRUE)
mean_lambda = 0.7*mean(eigenvalue)
dev.off()


