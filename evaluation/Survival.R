library(survminer)   # Load survminer first
library(survival)    # Then load survival
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

data <- read.csv(file = paste0(dirDataset, "matrix.csv"), row.names = 1)

metadata <- read.csv(file = paste0(dirDataset, "metadata.csv"), row.names = 1)

# Step 1: Split metadata based on the `source_name_ch1` column
list <- split(metadata$geo_accession, metadata$`source_name_ch1`)

# Step 2: Retrieve healthy controls (Non-tumor) and cases (Tumor)
control <- list$`Non-tumor`
case <- list$`Tumor`

dataC <- data[, case]     # Subset for case samples (Tumor)

print(metadata)
print(dataC)

###########GENE
gene <- "hsa-miR-100"

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

# Ensure that case_id and geo_accession are in the same format (character type)
df$case_id <- as.character(df$case_id)
metadata$geo_accession <- as.character(metadata$geo_accession)

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
