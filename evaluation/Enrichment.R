# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)
################################################
getEnrichmentPlot <- function(annotation,type,top_term,thr_pval,dirEnrich){

  annotation <- annotation[annotation$Adjusted.P.value < thr_pval, ]
  annotation <- annotation[order(annotation$Adjusted.P.value, decreasing = F),]

  annotation$Gene_count <- sapply(annotation$Genes, function(x){

    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)

  })

  annotation$Gene_ratio <- unlist(lapply(annotation$Overlap, function(x){

    total <- as.numeric(strsplit(x,"/")[[1]][2])
    count <- as.numeric(strsplit(x,"/")[[1]][1])

    Gene_ratio <- count/total

  } ))

  if(length(top_term) != 0 & top_term <= nrow(annotation)){
    annotation_top <- annotation[1:top_term,]
  }else{
    annotation_top <- annotation
  }

  # bar plot
  g1 <- ggplot(annotation_top,
               aes_string(x = "Gene_count", y = fct_reorder(annotation_top$Term, annotation_top$Gene_count), fill = "Adjusted.P.value" )) +
    geom_bar(stat="identity") +
    scale_fill_continuous(low="red", high="blue", name = "Adjusted.P.value",
                          guide=guide_colorbar(reverse=TRUE)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme_bw(base_size = 10) +
    ylab(NULL)

  pdf(paste0(dirEnrich,type,"_barplot.pdf"))
  print(g1)
  dev.off()

  # dot plot
  g2 <- ggplot(annotation_top,
               aes_string(x = "Gene_count", y = fct_reorder(annotation_top$Term, annotation_top$Gene_count))) +
    geom_point(aes(size = Gene_ratio, color = Adjusted.P.value )) +
    scale_colour_gradient(limits=c(min(annotation_top$Adjusted.P.value), max(annotation_top$Adjusted.P.value)), low = "red", high = "blue") +
    theme_bw(base_size = 10) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ylab(NULL)

  pdf(paste0(dirEnrich,type,"_dotplot.pdf"))
  print(g2)
  dev.off()

  write.table(annotation[,c("Term","Overlap","P.value","Adjusted.P.value","Gene_count","Gene_ratio","Genes")],
              paste0(dirEnrich,type,"_adj_pval_",thr_pval,".txt"), sep = "\t", quote = F, col.names = T, row.names = F)


}

################################################
setwd("~/")
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

dirEnrich <- paste0(dirDataset,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}
################################################
top_term <- 10
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirDataset,"DEG.txt")

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "DisGeNET")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
#input_list <- input_list$genes
input_list$genes <- sub("///.*", "", input_list$genes)

list <- split(input_list$genes, input_list$direction)

df_UP <- enrichr(list$UP, dbs)
df_DOWN <- enrichr(list$DOWN, dbs)

BP_UP <- df_UP$GO_Biological_Process_2021
MF_UP <- df_UP$GO_Molecular_Function_2021
DisGeNET_UP <- df_UP$DisGeNET
KEGG_UP <- df_UP$KEGG_2021_Human

BP_DOWN <- df_DOWN$GO_Biological_Process_2021
MF_DOWN <- df_DOWN$GO_Molecular_Function_2021
DisGeNET_DOWN <- df_DOWN$DisGeNET
KEGG_DOWN <- df_DOWN$KEGG_2021_Human

getEnrichmentPlot(BP_UP,"GO_BP_UP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(MF_UP,"GO_MF_UP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(KEGG_UP,"KEGG_UP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(DisGeNET_UP,"DisGeNET_UP",top_term,thr_pval,dirEnrich)

getEnrichmentPlot(BP_DOWN,"GO_BP_DOWN",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(MF_DOWN,"GO_MF_DOWN",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(KEGG_DOWN,"KEGG_DOWN",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(DisGeNET_DOWN,"DisGeNET_DOWN",top_term,thr_pval,dirEnrich)

