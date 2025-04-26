rm(list = ls())
options(stringsAsFactors = F)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")


library(stringr)
library(pheatmap) #heatmap

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
# now i have to import my data matrix
# path = "Results/GEO_HCC/" #se vogliamo sostituire con il dataset di TCGA il pat diventa vuoto
filename_in = paste0(dirDataset,"matrix.txt") # la matrice per TCGA è matrix_RNAseq_BRCA.txt

######################
#file in output

filename_DEG = paste0(dirDataset,"DEG.txt") # the most important file of differentials expressed genes
filename_list_normal= paste0(dirDataset, "normal.txt")
filename_list_case= paste0(dirDataset, "case.txt")
filename_DEG= paste0(dirDataset, "DEG.txt")
filename_var= paste0(dirDataset, "var.txt")
filename_matrix_DEG= paste0(dirDataset, "matrix_DEG.txt")
filename_heatmap= paste0(dirDataset, "heatmap.png")

######################
prc_IQR = 0.1
thr_fc = 1.2
thr_pval = 0.05
paired = FALSE
#######################
# STEP 1 : importing the data from Geoquery
#######################

readr::local_edition(1)
library("GEOquery")

setwd("~/")
##### Set Results folder ####
dirRes <- "Results/"
if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory ",dirRes," already exists"))
}
dataset <- "RA"
dirDataset <- paste0(dirRes,dataset,"/")
if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory ",dirDataset," already exists"))
}

###########################################
# STEP 1: Downloading data
###########################################
series = "GSE10072" #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93272
tmp <- getGEO(GEO = series) #show(tmp)
set <- tmp$GSE10072_series_matrix.txt.gz
pData <- phenoData(set)
metadata <- pData@data
aData <- assayData(set)
matrix <- data.frame(aData$exprs)
matrix<- 2^matrix
rm(pData,aData,tmp)

###########################################
# STEP 2: Preparing data
###########################################
annotation <- fData(set)
geneSymbol <- annotation$`Gene Symbol`

matrix <- matrix[annotation$ID,]
matrix <- aggregate(matrix,list(geneSymbol),"mean")
ind <- which(matrix$Group.1 == "")
if(length(ind) > 0) matrix <- matrix[-ind,]
rownames(matrix) <- matrix$Group.1
matrix <- matrix[,-1]
matrix[is.na(matrix)] <- 0
rm(set,annotation,ind,geneSymbol)

############################################
# STEP 3: Extracting case and control samples
###########################################
list <- split(metadata$geo_accession,
              metadata$`source_name_ch1`)
# We have to take healthy control and RA
control <- list$`Normal Lung Tissue`
case <- list$`Adenocarcinoma of the Lung` # for retrieving instead MCI data (case <- list$MCI)
data <- matrix[,c(control,case)]
dataN <- matrix[,control]
dataC <- matrix[,case]
rm(list,matrix)

############################################
#STEP 4: Export data
###########################################
write.table(data,
            paste0(dirDataset,"matrix.txt"), sep= "\t",
            col.names = NA, row.names = T, quote = F)
write.table(control,paste0(dirDataset,
                           "normal.txt"), sep= "\t", col.names = F,
            row.names = F, quote = F)
write.table(case,
            paste0(dirDataset,"case.txt"), sep= "\t",
            col.names = F, row.names = F, quote = F)
write.table(metadata,paste0(dirDataset,"metadata.txt"), sep= "\t", col.names = T,
            row.names = F, quote = F)
# Remove downloaded data
unlink(series, recursive = TRUE)

#################
#step 2: ANALYSIS
#################
# overall mean
genes <- rownames(data)

overall_mean = apply (data,  1, mean ) #we eliminate those one with overall 0
ind= which(overall_mean == 0)
if (length(ind) > 0){   #check to see if ind is not empty
  dataN = dataN[-ind,]
  dataC = dataC[-ind,]
  data = data[-ind,]
  genes = genes[-ind]
}

#step 2.1 : logarithmic trasformation
dataN = log2(dataN+1)
dataC = log2(dataC+1)
data = log2(data+1)

write.table(dataN,
            paste0(dirDataset,"dataN_Log.txt"), sep= "\t",
            col.names = NA, row.names = T, quote = F)


write.table(dataC,
            paste0(dirDataset,"dataC_Log.txt"), sep= "\t",
            col.names = NA, row.names = T, quote = F)

#step 2.2 : pre-processing
variation = apply (data, 1, IQR) # for each gene we apply the function IQR.
# IQR help us to measure the variation around the median of our genes

variat_data <- data.frame(Gene = genes, variation = variation)
write.table(variat_data, file = filename_var, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Print a message with the file path
cat("variation values saved to file: ", filename_var)


#IQR filtering
thr_prc = quantile(variation,prc_IQR) # we quantify a percentile (thr_prc with prc_iqr=0.1 is 0.1592515)
ind = which (variation <= thr_prc)
dataN = dataN[ -ind,]
dataC = dataC[ -ind,]
data= data[-ind,]
genes = genes[-ind]
rm(ind)

#devo scegliere

#iqr distribution
hist ( variation,
       main = "IQR Frequency distribution", breaks = 100 ,xlab= "IQR value",ylab="Frequency", col = "seagreen")
abline(v= thr_prc, lty = 2 , lwd = 4 , col="darkgoldenrod1")

# step 2.3 Log fc

#logFC > 0 is increasing its expression in case samples
#logfc < 0 is decreasing its expression in case samples
logFC = rowMeans(dataC) - rowMeans(dataN)
print(log2(thr_fc))
hist(logFC, main = " FC (logarithmic) frequency distribution",breaks = 100,xlab = " logFC",ylab = " frequency",col ="orange")
abline(v = c(-log2(thr_fc), log2(thr_fc)),lty =2 ,lwd = 4, col = "red")

ind = which(abs(logFC) < log2(thr_fc)) # we are eliminating all those genes that are satisfying this condition

if(length(ind)>0)
{
  dataN = dataN[-ind,]
  dataC = dataC[-ind,]
  data = data[-ind,]
  genes = genes[-ind]
  logFC = logFC [-ind]

}

rm(ind)

#now we have to assign a statistical significance to those changes

# step 2.4 p-value computation
N = ncol(dataN)
M = ncol(dataC)
#x is the first gene of the data matrix
pval = apply(data,1,function(x){
  t.test(x[1:N],x[(N+1):(M+N)],paired=paired)$p.value
})

#adjustment p-value
pval_adj = p.adjust(pval,method="fdr")

#p-value filtering
ind = which(pval_adj > thr_pval)
if(length(ind)>0)
{
  dataN = dataN[-ind,]
  dataC = dataC[-ind,]
  data = data[-ind,]
  genes = genes [-ind]
  logFC = logFC[-ind]
  pval = pval [-ind]
  pval_adj = pval_adj[-ind]
}
rm(ind)

# Volcano plot
# Define the color and shape for each point
color <- ifelse(pval_adj > 0.05 | abs(logFC) <= log2(thr_fc), "black", ifelse(logFC > 0 & pval_adj < 0.05, "red", "green"))

# Create the volcano plot
plot(logFC, -log10(pval_adj), pch = 16, col = color, cex = 0.7, main = "Volcano Plot", xlab = "logFC", ylab = "-log10(Adjusted P-value)")

# Add horizontal line at p-value threshold
abline(h = -log10(0.05), col = "blue", lty = 2)

# Add vertical lines at logFC threshold
abline(v = c(-log2(thr_fc), log2(thr_fc)), col = "blue", lty = 2)


##############
#STEP 3 : exporting results
##############
results <- data.frame(genes = genes, logFC = logFC, pval = pval, pval_adj = pval_adj)
direction <- ifelse(logFC > 0, "UP", "DOWN")
results <- cbind(results, direction)

#we are splitting in 2 pieces
#DEG = data.frame (str_split_fixed(genes,"\\|",2))  #non serve
# colnames(DEG)=c("geneSymbol","ensembl_id")


write.table(results,file = filename_DEG, row.names = F, sep = "\t", quote =F)
write.table(data, file = filename_matrix_DEG,row.names = T,col.names = NA, sep="\t", quote = F)
#write.table(pzC_com, filename_list_case,sep = "\t ", col.names = F, row.names = F,quote = F)

###########################
# STEP 4 : PLOTS
###########################

plot(logFC, -log10(pval),
     main = "Volcano Plot",
     xlim = c(-4, 4),
     ylim = c(0, 70),
     xlab = "log2 fold change",
     ylab = "-log10 p-value"
)

abline(h= - log10(thr_pval), lty = 2, lwd= 4, col= "red")
abline(v = c(-log2(thr_fc), log2(thr_fc)),lty= 2,lwd = 4, col="grey")

#Pie chart
count <- table(results$direction)
pie(count,
    labels= paste0(names(count), " ", round(100 * count/sum(count),2), "%"),
    col = c("blue","gold"))

# box plot
ind = which.max(logFC) #most up-regulated
gene_id = results$genes[ind]

boxplot(as.numeric(dataN[ind,]), as.numeric(dataC[ind,]),
        main=paste(gene_id,", ", "adjusted p-value= ",format(pval_adj[ind],digits =2)),
        notch = T,
        ylab="Gene expression value",
        xlab= "Condition",
        col= c("green","violet"),
        names = c("Normal", "Case"),
        pars= list(boxwex = 0.3, staplewex = 0.6))

rm(ind)

ind = which.min(logFC) #most down-regulated
gene_id = results$genes[ind]

boxplot(as.numeric(dataN[ind,]), as.numeric(dataC[ind,]),
        main=paste(gene_id,", ", "adjusted p-value= ",format(pval_adj[ind],digits =2)),
        notch = T,
        ylab="Gene expression value",
        xlab= "Condition",
        col= c("green","violet"),
        names = c("Normal", "Case"),
        pars= list(boxwex = 0.3, staplewex = 0.6))

rm(ind)


df1 = data.frame(logFC = logFC, pval = -log10(pval_adj))
condition1 <- (pval_adj <= thr_pval) & (logFC > log2(thr_fc))
condition2 <- (pval_adj <= thr_pval) & (logFC < -log2(thr_fc))
df1$legend <-ifelse(condition1,"up",ifelse(condition2,"down","delete"))

p <- ggplot(df1, aes(x = logFC, y = pval, color = legend)) +
  geom_point() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("delete" = "darkgrey", "up" = "red", "down" = "green")) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "Volcano plot", x = "Log2 Fold change (FC)", y = "Log10( adjusted p value )") +
  geom_hline(yintercept = -log10(0.01), linetype = 2, color = "black") +
  geom_vline(xintercept = log2(thr_fc), linetype = 2, color = "black") +
  geom_vline(xintercept = -log2(thr_fc), linetype = 2, color = "black")

print(p)

###########################
#HEAT:MWP
#
annotation <- data.frame(condition = metadata$`characteristics_ch1.1`)
rownames(annotation) <- metadata$geo_accession

vect_color <- c("green","violet")
names(vect_color) <- unique(annotation$condition)

annotation_colors <- list(condition = vect_color)
# Stiamo plottando i dati diespressione

pheatmap(data, scale = "row",
         border_color = NA,
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_rows = "correlation", #clustering con correlation (similarity)
         clustering_distance_cols = "correlation",
         clustering_method = "average", # dato l'insieme sono vicini con il valor medio
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         color = colorRampPalette(colors = c("blue", "blue3", "black", "yellow3", "yellow"))(100),
         show_rownames = F,
         show_colnames = F,
         cutree_cols = 2,
         cutree_rows = 2,
         width = 10, height = 10,
         filename = filename_heatmap)

