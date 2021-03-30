###########################################################
#   Script to run scrna-seq preprocessing workflow        #
#               D.Croucher (L.Richards template)          #
#                      Sept 2020                          #
###########################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Collect sequencing metrics outputted by cellranger count
### 2) Perform cell barcode calling using dropbead
### 3) Read data and merge into single seurat object 
### 4) Add metadata
### 5) Add cell-level quality control metrics 
### 6) Add QC flag
### 7) Output plots
### 8) Save data

###########################################################

###########################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 72:00:00
## #SBATCH --mem=60G
## #SBATCH -p himem
## #SBATCH -c 15
## #SBATCH -N 1
##
## module load R/3.6.1
##
## Rscript /cluster/projects/pughlab/projects/scVKMYC/scripts/preprocessing.R --samples /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/vkmyc_cohort_sampleInput.txt \
##    --cellrangerOutputFolder /cluster/projects/pughlab/projects/scVKMYC_mPC/cr_outs/ \
##    --bamTagHistogramOutputFolder /cluster/projects/pughlab/projects/scVKMYC_mPC/BAMTagHistogram/ \
##    --outputDir /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/seurat/ \
##    --rawDGEinputFolder /cluster/projects/pughlab/projects/scVKMYC_mPC/cr_outs/ \
##    --sampleShortNames /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/vkmyc_cohort_sampleshortInput.txt \
##    --fileName scVKMYC_Cohort_preprocess \
##    --sampleMetaData /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/vkmyc_cohort_sampleMetaData.csv \
##    --outputPlots TRUE 
##
###########################################################

### start clocks
StopWatchStart <- list()
StopWatchEnd   <- list()
StopWatchStart$Overall  <- Sys.time()
StopWatchEnd$Overall  <- Sys.time()

#########################################
# 1) Load packages
#########################################

print("")
print("********************")
print("Load packages")
print(Sys.time())
print("********************")
print("")
StopWatchStart$LoadPackages  <- Sys.time()

suppressMessages(library(scran))
suppressMessages(library(Matrix))
suppressMessages(library(scales))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(DropletUtils))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(sctransform))
suppressMessages(library(future))
suppressMessages(library(KneeArrower))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(grid))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
suppressMessages(library(gplots))
suppressMessages(library(heatmap.plus))
suppressMessages(library(readxl))
suppressMessages(library(Rmisc))
suppressMessages(library(stringr))
suppressMessages(library(devtools))
suppressMessages(library(plyr))
suppressMessages(library(dropbead))


StopWatchEnd$LoadPackages  <- Sys.time()


#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--samples",
                                type = "character",
                                default = NULL,
                                help = "path to txt file with column SAMPLEID; each sample on own line",
                                metavar= "character"
                               ),
                     make_option("--fileName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files, type NONE is you want to use sampleID",
                                metavar= "character"
                               ),
                     make_option("--cellrangerOutputFolder",
                                type = "character",
                                help = "path to dir with cr count outs to read in",
                                default = NULL,
                                metavar= "character"
                               ),
                     make_option("--bamTagHistogramOutputFolder",
                                type = "character",
                                help = "path to dir with bamtaghistogram outs to read in",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--outputDir",
                                type = "character",
                                help = "path to dir where all pipeline outputs will go",
                                default = NULL,
                                metavar= "character"
                               ),
                     make_option("--rawDGEinputFolder",
                                type = "character",
                                help = "same as cellrangerOutputFolder",
                                default = NULL,
                                metavar= "character"
                               ),
                     make_option("--outputPlots",
                                type = "logical",
                                help = "TRUE/FALSE to output plots",
                                default = NULL,
                                metavar= "logical"
                               ),                     
                    make_option("--sampleMetaData",
                                type = "character",
                                help = "path to dir where sample level metadata file is located",
                                default = NULL,
                                metavar= "character"
                               ),
                     make_option("--sampleShortNames",
                                type = "character",
                                help = "path to txt file with column SAMPLEID; each sample on own line",
                                default = NULL,
                                metavar= "character"
                               )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

samples <- opt$samples
fileName <- opt$fileName
cellrangerOutputFolder <- opt$cellrangerOutputFolder
bamTagHistogramOutputFolder <- opt$bamTagHistogramOutputFolder
rawDGEinputFolder <- opt$rawDGEinputFolder
sampleShortNames <- opt$sampleShortNames
outputDir <- opt$outputDir
outputPlots <- opt$outputPlots
sampleMetaData <- opt$sampleMetaData

#########################################
# Set up parallization
#########################################

print("")
print("********************")
print("Setup Parallelization with future")
print(Sys.time())
print("********************")
print("")
StopWatchStart$Parallelize  <- Sys.time()

CoresAvailable <- as.numeric(availableCores()[[1]])
print(paste0(CoresAvailable, " Cores Available"))
plan("multicore", workers = CoresAvailable)
options(future.globals.maxSize = 10 * 1024 ^ 3)
#plan()

StopWatchEnd$Parallelize  <- Sys.time()


#########################################
# Pipeline Variables
#########################################

print("")
print("********************")
print("Pipeline Variables")
print(Sys.time())
print("********************")
print("")
StopWatchStart$SetUpVaribles  <- Sys.time()

sampleNames <- c(as.character(read.table(samples, header = T)$SAMPLEID))
if(length(sampleNames) == 1){
    print(paste0("Sample IDs: ", sampleNames))
} else if(length(sampleNames > 1)){
    print(paste0("Sample IDs: ", paste(sampleNames, collapse = ", ")))
}

sampleShortNames <- c(as.character(read.table(sampleShortNames, header = T)$SAMPLEID))
if(length(sampleShortNames) == 1){
    print(paste0("Sample IDs: ", sampleShortNames))
} else if(length(sampleShortNames > 1)){
    print(paste0("Sample IDs: ", paste(sampleShortNames, collapse = ", ")))
}
  
dataDir <- c(paste(cellrangerOutputFolder, 
                   sampleNames, 
                   "/raw_gene_bc_matrices/testing/", 
                   sep = ""))
  
  
if(fileName == "NONE"){assign("fileName", sampleNames)}
print(paste0("File prefix: ", fileName))


StopWatchEnd$SetUpVaribles  <- Sys.time()

setwd(outputDir)
dir.create(fileName)
setwd(fileName)
dir.create("data")
if(outputPlots == TRUE){dir.create("figures")}


#########################################
# Collect sequencing metrics outputted by cellranger count
#########################################

print("")
print("*********************")
print("Collect sequencing metrics")
print(Sys.time())
print("********************")
print("")
StopWatchStart$SeqMetrics <- Sys.time()


seqMetrics <- data.frame(rep("", 20)) # 20 = number of metrics in summary.csv doc

for (i in seq_along(sampleNames)) {
  
  seqMetrics[,i] <- t(read.csv(paste(cellrangerOutputFolder, 
                                     sampleNames[i], 
                                     "/metrics_summary.csv", sep = "")))

}

colnames(seqMetrics) <- sampleNames
rownames(seqMetrics) <- rownames(t(read.csv(paste(cellrangerOutputFolder, 
                                     sampleNames[i], 
                                     "/metrics_summary.csv", sep = ""))))

seqMetrics.file <- paste0("./data/", fileName, "_seqMetrics.csv")
write.csv(seqMetrics, file = seqMetrics.file)


StopWatchEnd$SeqMetrics <- Sys.time()




#########################################
# Perform cell barcode calling using dropbead
#########################################

print("")
print("*********************")
print("Cell barcode calling")
print(Sys.time())
print("********************")
print("")
StopWatchStart$CellBarcodeCalling <- Sys.time()

#Read in readcount files

readCounts <- list()

for (i in seq_along(sampleNames)) {
  
  readCounts[[i]] <- read.table(paste(bamTagHistogramOutputFolder,
                                     sampleNames[i],
                                     ".bam.readcounts.txt.gz", 
                                     sep = ""), 
                               header = F, 
                               stringsAsFactors = FALSE)
  readCounts[[i]][[1]] <- as.numeric(readCounts[[i]][[1]])
}

names(readCounts) <- sampleNames


#Perform barcode calling

barcodeTotal <- data.frame(rep("", 1))

for (i in seq_along(sampleNames)) {

  barcodeTotal[,i] <- dropbead::estimateCellNumber(readCounts[[i]][[1]], 
                                          max.cells = 50000)

  }

colnames(barcodeTotal) <- sampleNames


##Collect cell barcodes and re-name to match what will be incomming dge barcode names (colnames)

barcodeCalls <- list()

for (i in seq_along(sampleNames)) {
  
  barcodesToUse <- barcodeTotal[1,i]
  tmpPaste <- rep(sampleShortNames[i], barcodesToUse) 
  tmpBarcodes <- (readCounts[[i]][[2]])[1:barcodesToUse] 
  barcodeCalls[[i]] <- paste(tmpPaste, tmpBarcodes, sep = "_")
}

barcodeCalls <- lapply(X = barcodeCalls, function(x) gsub("-1", "", x)) #remove -1
names(barcodeCalls) <- sampleNames


StopWatchEnd$CellBarcodeCalling <- Sys.time()

              
                       
#########################################
# Read data and merge into single seurat object 
#########################################

print("")
print("*********************")
print("Read data and merge")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadDataMerge <- Sys.time()

#Load and merge DGEs into a single seurat object

for (i in seq_along(sampleShortNames)) {
  
  seuratTmp <- Read10X(dataDir[i], strip.suffix = TRUE) #returns the full raw_dge matrix with cellnames updated to ND1_XXXXXXXXX
  colnames(seuratTmp) <- paste(sampleShortNames[i], colnames(seuratTmp), sep="_") 
  
  seuratTmp <- seuratTmp[, colnames(seuratTmp) %in% barcodeCalls[[i]]] #subset to only include cells called from dropbead
  
  
  seuratTmp <- CreateSeuratObject(counts = seuratTmp)
  
  
    if (i == as.integer(1)){

      seurat <- seuratTmp
    
    }
  
    if (i > as.integer(1)){
      
      seurat <- merge(x=seurat, y=seuratTmp)
    }

}

print("Total Dataset Size")
print(dim(seurat@assays$RNA@data)) ## print out final size of object                       

StopWatchEnd$ReadDataMerge <- Sys.time()
                       
                       
                       
                       
#########################################
# Add metadata
#########################################

print("")
print("*********************")
print("Add metadata")
print(Sys.time())
print("********************")
print("")
StopWatchStart$AddMetaData <- Sys.time()


sample.meta.data <- read.csv(sampleMetaData)

samples <- unique(seurat@meta.data$orig.ident) #vector of sample names which is added to "orig.ident" in previous steps
meta.name <- colnames(sample.meta.data) #name of sample-level metadata (columns)

#set sample.meta.data row order to match samples

sample.meta.data$Sample_ID <- factor(sample.meta.data$Sample_ID, levels = samples)
sample.meta.data <- sample.meta.data[order(sample.meta.data$Sample_ID), ]

meta <- seurat@meta.data

#lift sample-level metadata over to each cell (row in meta.data)
for (i in seq_along(meta.name)) {
  
  for (j in seq_along(samples)) {
    
    meta[meta$orig.ident %in% samples[j], meta.name[i]] <- sample.meta.data[j,i]
    
  }
  
  #set factor levels based on order or sample-level metadata entries 
  meta[,meta.name[i]] <- factor(meta[,meta.name[i]], levels = unique(sample.meta.data[,i]))
  
}


seurat <- AddMetaData(seurat, meta)

StopWatchEnd$AddMetaData <- Sys.time()
                       
                       
                       
#########################################
# Add cell-level quality control metrics 
#########################################

print("")
print("*********************")
print("Add QC metrics")
print(Sys.time())
print("********************")
print("")
StopWatchStart$AddQC <- Sys.time()


###----ReadCounts

for (i in seq_along(sampleShortNames)) {
  readCounts[[i]]$V2 <- paste(sampleShortNames[i], 
                              readCounts[[i]]$V2, 
                              sep = "_")
  readCounts[[i]]$V2 <- gsub("-1", "", readCounts[[i]]$V2)
}

#subset readCounts to be only as long as barcode Total

for (i in seq_along(sampleShortNames)) {
  readCounts[[i]] <- readCounts[[i]][c(1:barcodeTotal[, i]), ]
}

#rbind all readCounts into one readCountMeta df

for (i in seq_along(sampleShortNames)) {
  
  if (i == 1) {
    readCountMeta <- as.data.frame(readCounts[[i]], 
                                   row.names = readCounts[[i]]$V2) 
    readCountMetaTmp <- readCountMeta
  }
  
  if (i > 1) {
    readCountMetaTmp <- as.data.frame(readCounts[[i]],
                                      row.names = readCounts[[i]]$V2) 
    readCountMeta <- rbind(readCountMeta, readCountMetaTmp)
  }
}

rm(readCountMetaTmp)
readCountMeta$V2 <- NULL
colnames(readCountMeta) <- "readCounts"

#Add in readCounts metadata
seurat <- AddMetaData(seurat, readCountMeta, col.name = "readCounts")


#Add mito metadata
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

#Fix nUMI/nGene
seurat[["nGene"]] <- seurat[["nFeature_RNA"]]
seurat[["nUMI"]] <- seurat[["nCount_RNA"]]              


StopWatchEnd$AddQC <- Sys.time()

                       
                       
#########################################
# Add QC flag
#########################################

print("")
print("*********************")
print("Add QC flag")
print(Sys.time())
print("********************")
print("")
StopWatchStart$QCFlag <- Sys.time()
                       
                       
###Flag outliers of percent.mt 

mitoQC <- WhichCells(seurat, expression = percent.mt > 15)
seurat@meta.data[rownames(seurat@meta.data) %in% mitoQC, "mitoQC"] <- "FAIL"
seurat@meta.data[is.na(seurat@meta.data$mitoQC), "mitoQC"] <- "PASS"


###Flag outliers of nGene

geneQC <- WhichCells(seurat, expression = nGene < 500)
seurat@meta.data[rownames(seurat@meta.data) %in% geneQC, "geneQC"] <- "FAIL"
seurat@meta.data[is.na(seurat@meta.data$geneQC), "geneQC"] <- "PASS"


###Flag outliers of nUMI

umiQC <- WhichCells(seurat, expression = nUMI < 1000)
seurat@meta.data[rownames(seurat@meta.data) %in% umiQC, "umiQC"] <- "FAIL"
seurat@meta.data[is.na(seurat@meta.data$umiQC), "umiQC"] <- "PASS"

###Create Cell-level QC flag (Cell_QC)

a <- paste0(seurat@meta.data$nGene_QC,
            seurat@meta.data$nUMI_QC,
            seurat@meta.data$mitoQC)
seurat@meta.data$Cell_QC <- ifelse(grepl("FAIL",a),
                                       "FAIL",
                                       "PASS"
                                      )                       

StopWatchEnd$QCFlag <- Sys.time()
                       
                       
#########################################
# Output plots
#########################################

if(outputPlots == TRUE){

    #############
    print("")
    print("********************")
    print("Output Plots")
    print(Sys.time())
    print("********************")
    print("")
    StopWatchStart$Plotting <- Sys.time()



    meta <- data.frame(seurat@meta.data)

    #######################

    print("1/6....QC plots")

    QC1 <- ggplot(meta, aes(x=nGene, y=nUMI)) +
          geom_point(size = 0.6, alpha = 0.7) +
          theme(legend.position="none") + theme_classic()
    QC1 <- ggMarginal(QC1, type="histogram")

    QC2 <- ggplot(meta, aes(x=nGene, y=percent.mt)) +
          geom_point(size = 0.6, alpha = 0.7) +
          theme(legend.position="none") + theme_classic()
    QC2 <- ggMarginal(QC2, type="histogram")


    qc.plots <- paste0("./figures/",fileName, "_QC.pdf")
    print(qc.plots)
    pdf(qc.plots, width = 10, height = 5)
    print(ggarrange(QC1, QC2, ncol = 2, nrow = 1))
    dev.off()

    
    StopWatchEnd$Plotting <- Sys.time()

    }







#########################################
# Save Data
#########################################

print("")
print("********************")
print("Save Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$SaveData <- Sys.time()

print("Raw Counts....")
raw.count.dir <- paste0("./data/", fileName, "_rawCounts")
raw <- Matrix(seurat@assays$RNA@counts, sparse = TRUE)
write10xCounts(path = raw.count.dir,
               x = raw
               )

print("Normalized Counts....")
norm.count.dir <- paste0("./data/", fileName, "_normCounts")
norm.count <- Matrix(seurat@assays$RNA@data, sparse = TRUE)
write10xCounts(path = norm.count.dir,
               x = norm.count
               )

print("Meta Data....")
meta <- data.frame(seurat@meta.data)
meta.file <- paste0("./data/", fileName, "_metaData.csv")
write.csv(meta, file = meta.file)
meta.file.2 <- paste0("./data/", fileName, "_metaData.rds")
saveRDS(meta, file = meta.file.2)



print("Saving Seurat Object....")
seurat.file <- paste0("./data/", fileName, "_seurat.rds")
saveRDS(seurat, file = seurat.file)

StopWatchEnd$SaveData <- Sys.time()
StopWatchEnd$Overall <- Sys.time()


#########################################
# Stop watches
#########################################

print("")
print("********************")
print("Run Time")
print(Sys.time())
print("********************")
print("")

clocks <- mapply('-', StopWatchEnd, StopWatchStart, SIMPLIFY = FALSE)
print(clocks)
#bb <- as.numeric(bb, units = "mins")
#lapply(bb, function(x){as.numeric(x, units = "mins")})


#########################################
# Session Information
#########################################

print("")
print("********************")
print("Session Info")
print(Sys.time())
print("********************")
print("")

print(sessionInfo())

print("***************************")
print("*******END OF SCRIPT*******")
print("***************************")
