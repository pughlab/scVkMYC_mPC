###########################################################
#   Script to predict doublets using DoubletFinder        #
#           D.Croucher (L.Richards template)              #
#                      Sept 2020                          #
###########################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load QC'ed data
### 2) Run doubletFinder pipeline
### 3) Flag doublets
### 4) Output plots
### 5) Save results


###########################################################

###########################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 10:00:00
## #SBATCH --mem=60G
## #SBATCH -p himem
## #SBATCH -c 15
## #SBATCH -N 1
## #SBATCH --account=pughlab
##
## module load R/3.6.1
##
## Rscript /cluster/projects/pughlab/projects/scVKMYC/scripts/doubletFinder.R \
## --seuratObj /cluster/projects/pughlab/projects/scVKMYC/analysis/test/seurat/scVKMYC_test/data/scVKMYC_test_seurat.rds \
## --outName scVKMYC_test_doubletFlagged \
## --outputDir /cluster/projects/pughlab/projects/scVKMYC/analysis/test/doubletFinder \
## --knownMultRate /cluster/projects/pughlab/projects/scVKMYC/analysis/KnownMultipletRate.csv \
## --numVarGenes 3000
##
###########################################################

### start clocks
StopWatchStart <- list()
StopWatchEnd   <- list()
StopWatchStart$Overall  <- Sys.time()
StopWatchEnd$Overall  <- Sys.time()

#########################################
# Load packages
#########################################

print("")
print("********************")
print("Load packages")
print(Sys.time())
print("********************")
print("")
StopWatchStart$LoadPackages  <- Sys.time()

suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(heatmap.plus))
suppressMessages(library(readxl))
suppressMessages(library(Rmisc))
suppressMessages(library(stringr))
suppressMessages(library(Hmisc))
suppressMessages(library(devtools))
suppressMessages(library(plyr))
suppressMessages(library(DoubletFinder))
suppressMessages(library(dropbead))
suppressMessages(library(optparse))



#########################################
# Functions
#########################################

#-------------Sample-specific doublet rate

# Zheng GXY et al. 2017. A cell titration experiment across six different cell loads showed a linear relationship between 
# the multiplet rate and the number of recovered cells ranging from 1,200 to 9,500 (Supplementary Fig. 1a).-- based on 
# version 1 chemistry -- but the point is its linear so assuming it can be used for computing doublet rate of other chemistries
# https://assets.ctfassets.net/an68im79xiti/RT8DYoZzhDJRBMrJCmVxl/6a0ed8015d89bf9602128a4c9f8962c8/CG00052_SingleCell3_ReagentKitv2UserGuide_RevF.pdf
# MultipletRate,RecoveredCells
# 0.4,500
# 0.8,1000
# 1.6,2000
# 2.3,3000
# 3.1,4000
# 3.9,5000
# 4.6,6000
# 5.4,7000
# 6.1,8000
# 6.9,9000
# 7.6,10000

predictDoubletRate <- function(object) {
  n.cells <- ncol(object)
  mod=lm(MultipletRate ~ RecoveredCells, data=known.multiple.rate)
  sample.multiple.rate = coef(mod)[1] + n.cells*coef(mod)[2]
}



#-------------Modified function for predicting doublets (based on doubletFinder V2)

doubletFinderMod <- function (seu, expected.doublets = 0, proportion.artificial = 0.25, 
                              proportion.NN = 0.01) 
{
  if (expected.doublets == 0) {
    stop("Need to set number of expected doublets...")
  }
  
  print("Creating artificial doublets...")
  
  data <- GetAssayData(object = seu, slot = "counts")[, colnames(seu)]
  
  real.cells <- colnames(seu)
  
  n_real.cells <- length(real.cells)
  
  n_doublets <- round(n_real.cells/(1 - proportion.artificial) - 
                        n_real.cells)
  
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  
  data_wdoublets <- cbind(data, doublets)
  
  print("Creating Seurat object...")
    seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)
  
  print("Normalizing Seurat object...")
  seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets)
  
  print("Finding variable genes...")
  seu_wdoublets <- Seurat::FindVariableFeatures(seu_wdoublets, selection.method = "vst", nfeatures = numVarGenes)
  
  print("Scaling data...")
  seu_wdoublets <- Seurat::ScaleData(seu_wdoublets, features = VariableFeatures(object = seu_wdoublets))
  
  print("Running PCA...")
  seu_wdoublets <- Seurat::RunPCA(seu_wdoublets, features = VariableFeatures(object = seu_wdoublets, npcs = 50, verbose = FALSE))
  
  cell.names <- colnames(seu_wdoublets)
  nCells <- length(cell.names)
  
  print("Calculating PC distance matrix...")
  PCs <- seu@commands$RunTSNE@params$dims
  
  if (length(PCs) == 0) {
    stop("Need to run tSNE on original Seurat object...")
  }
  
  pca.coord <- Embeddings(object = seu_wdoublets, reduction = "pca")[, PCs]

  rm(seu_wdoublets)
  gc()
  
  dist.mat <- as.matrix(dist(pca.coord))
  dist.mat <- dist.mat[, -grep("X", colnames(dist.mat))]
  
  pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
  
  rownames(pANN) <- real.cells
  colnames(pANN) <- "pANN"
  
  k <- round(nCells * proportion.NN)
  for (i in 1:n_real.cells) {
    neighbors <- order(dist.mat[, i])
    neighbors <- neighbors[2:(k + 1)]
    neighbor.names <- rownames(dist.mat)[neighbors]
    pANN[i, 1] <- length(grep("X", neighbor.names))/k
  }
  
  seu$pANN <- pANN
  predictions <- as.data.frame(rep("Singlet", n_real.cells), 
                               ncol = 1, stringsAsFactors = FALSE)
  
  rownames(predictions) <- real.cells
  
  doublet.predictions <- rownames(seu@meta.data)[order(seu@meta.data$pANN, 
                                                       decreasing = TRUE)]
  
  doublet.predictions <- doublet.predictions[1:expected.doublets]
  predictions[doublet.predictions, ] <- "Doublet"
  colnames(predictions) <- "pANNPredictions"
  
  seu$pANNPredictions <- predictions
  return(seu)

}


StopWatchEnd$LoadPackages  <- Sys.time()


#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--seuratObj",
                                type = "character",
                                default = NULL,
                                help = "path to seurat object containing SingleR metadata annotations",
                                metavar= "character"
                               ),
                     make_option("--outName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files",
                                metavar= "character"
                               ),
                     make_option("--outputDir",
                                type = "character",
                                help = "path to dir where results will be outputted",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--knownMultRate",
                                type = "character",
                                help = "path to file containing known multiple rates",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--numVarGenes",
                                type = "integer",
                                default = NULL,
                                help = "number of varibale genes to identify",
                                metavar= "integer"
                               )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

seuratObj <- opt$seuratObj
outName <- opt$outName
outputDir <- opt$outputDir
knownMultRate <- opt$knownMultRate
numVarGenes <- opt$numVarGenes


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

print(paste0("File prefix: ", outName))
print(paste0("Input path: ", seuratObj))
print(paste0("Path to file with known multiplet rates: ", knownMultRate))
print(paste0("Number Variable Genes: ", numVarGenes))

StopWatchEnd$SetUpVaribles  <- Sys.time()

setwd(outputDir)
dir.create(outName)
setwd(outName)


#########################################
# Load data
#########################################

print("")
print("*********************")
print("Loading Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadData <- Sys.time()

# Seurat object containing SingleR predictions
seurat <- readRDS(seuratObj)

# File containing known multiplet rates
known.multiple.rate <- read.csv(knownMultRate)

print("")
print("*********************")
print("Known Multiplet Rate")
print("********************")
print("")

print(known.multiple.rate)

StopWatchEnd$ReadData <- Sys.time()



#########################################
# Run doubletFinder pipeline (on each sample)
#########################################

# https://github.com/chris-mcginnis-ucsf/DoubletFinder

sampleShortNames <- levels(seurat@meta.data$Sample_ID)

for (i in seq_along(sampleShortNames)) {
 
  tmp <- SubsetData(seurat, 
                    cells = grep(sampleShortNames[i], 
                                 colnames(seurat), 
                                 value = T))
  
  #based on 10x estimated doublets and linear relationship b/t cell number and doublet rate
  sample.spec.exp.doublets <- predictDoubletRate(tmp)/100
  
  tmp <- doubletFinderMod(tmp,
                          expected.doublets = sample.spec.exp.doublets*ncol(tmp), 
                          proportion.artificial = 0.25, 
                          proportion.NN = 0.01) 
  
  tmpMeta <- data.frame(tmp@meta.data[,"pANNPredictions"],
                        row.names = rownames(tmp@meta.data))
  
  if (i == 1) {
  doubletDF <- tmpMeta
  
  }else{
    doubletDF <- rbind(doubletDF, tmpMeta)
  }
  
}

colnames(doubletDF) <- "pANNPredictions"

seurat <- AddMetaData(seurat, doubletDF, 
                      "pANNPredictions")

print("")
print("*********************")
print("Summary of DoubletFinder Results")
print("********************")
print("")

print("Total mutiplets detected")
doubletCells <- rownames(seurat@meta.data)[seurat@meta.data$pANNPredictions == "Doublet"]
length(doubletCells)

print("Sample-Level Summary")
table(seurat@meta.data[, "pANNPredictions"], seurat@meta.data[, "Sample_ID"])


#########################################
# Flag doublets
#########################################

###Flag outliers of nUMI

doubletQC <- doubletCells
seurat@meta.data[rownames(seurat@meta.data) %in% doubletQC, "doubletQC"] <- "FAIL"
seurat@meta.data[is.na(seurat@meta.data$doubletQC), "doubletQC"] <- "PASS"

###Create Cell-level QC flag (Cell_QC)

a <- paste0(seurat@meta.data$nGene_QC,
            seurat@meta.data$nUMI_QC,
            seurat@meta.data$mitoQC,
            seurat@meta.data$doubletQC)
seurat@meta.data$Cell_QC <- ifelse(grepl("FAIL",a),
                                       "FAIL",
                                       "PASS"
                                      )  

#########################################
# Output plots
#########################################

# --------- Flag doublets
  
doubletAssigned <- doubletCells

seurat@meta.data[rownames(seurat@meta.data) %in% doubletAssigned, "doubletAssignment"] <- "doublet"
seurat@meta.data[is.na(seurat@meta.data$doubletAssignment), "doubletAssignment"] <- "singlet"

singletCells <- rownames(seurat@meta.data)[seurat@meta.data$doubletAssignment == "singlet"]

tmpseurat <- subset(seurat, cells = singletCells)



# --------- tsne of doublets

meta <- seurat@meta.data

tsne <- ggplot(meta, aes(tSNE_1, tSNE_2))
tsne <- tsne + geom_point(aes(colour=factor(doubletAssignment)), alpha=0.7, size=0.5) + 
               labs(x = "", y = "") +
               scale_colour_brewer(palette = "Set1") + 
               theme_classic() +
               theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
               guides(color = guide_legend(override.aes = list(size=3))) +
               ggtitle(paste0(outName, "TSNE - Predicted Multiplets", sep = " "))

umap <- ggplot(meta, aes(UMAP_1, UMAP_2))
umap <- umap + geom_point(aes(colour=factor(doubletAssignment)), alpha=0.7, size=0.5) + 
               labs(x = "", y = "") +
               scale_colour_brewer(palette = "Set1") + 
               theme_classic() +
               theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
               guides(color = guide_legend(override.aes = list(size=3))) +
               ggtitle(paste0(outName, "UMAP - Predicted Multiplets", sep = " "))

tsne.file <- paste0(outName, "_tSNE.pdf")
pdf(tsne.file, width = 10, height = 10)
print(tsne)
dev.off()

umap.file <- paste0(outName, "_UMAP.pdf")
pdf(umap.file, width = 10, height = 10)
print(umap)
dev.off()

#########################################
# Save results
#########################################

print("Saving Seurat Object....")
seurat.file <- paste0(outName, "_seurat.rds")
saveRDS(seurat, file = seurat.file)

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
