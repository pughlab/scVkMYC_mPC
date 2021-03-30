########################################################################
#   Script to run harmony plus multi-resolution clustering pipeline    #
#     L.Richards (multi-res clustering) (D.Croucher mods/harmony)      #
#                               April 2020                             #
########################################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads lineage-subsetted scRNA-seq data (seurat objects)
### 2) Re-normalize (LogNormalize, scran, scTransform) and re-scale expression matrix
### 3) Identifies variable genes (vst)
### 4) Cell cyle scoring
### 5) Run PCA and determine number of significant PCs (can auto-detect scree)
### 6) Run Harmony batch correction
### 7) Non-linear dimensionality reduction (tSNE, UMAP)
### 8) Cluster cells over range of 6 resolutions
### 9) Runs differential gene expression (DGE) over all clustering resolutions
### 10) Calculates silhouette width for clusters
### 11) Selects optimal clustering resolution based on DE-gene cutoff and max silhouette width
### 12) Saves data
### 13) Output plots
###########################################################

###########################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 72:00:00
## #SBATCH --mem=300G
## #SBATCH -p superhimem
## #SBATCH -c 50
## #SBATCH -N 1
## #SBATCH --account=pughlab
##
## module load R/3.6.1
##
## Rscript /cluster/projects/pughlab/projects/scVKMYC/scripts/harmony_multires_clustering_pipeline.R \
##     --seuratObj /cluster/projects/pughlab/projects/scVKMYC/analysis/test/lineages/scVKMYC_test_lineages/B_PC/scVKMYC_test_lineages_B_PC_seurat.rds  \
##     --fileName scVKMYC_test_multires_harmony \
##     --outputDir /cluster/projects/pughlab/projects/scVKMYC/analysis/lineages/scVKMYC_test_lineages/B_PC/ \
##     --gene.filtering 0.001 \
##     --normalizationMethod LogNormalize \
##     --varsRegress /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/varsRegress.txt \
##     --numVarGenes 3000 \
##     --pcaScalingGenes var \
##     --computePC 75 \
##     --numPC 0 \
##     --addPC 0 \
##     --outputPlots TRUE \
##     --minResolution 1 \
##     --maxResolution 2 \
##     --deGeneCutoff 15 \
##     --fdrCutoff 0.05 \
##     --optimalCluster TRUE \
##     --markerGenes /cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/CellTypeMarkerGenes_B_PC_Cells.csv \
##     --species mouse
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
suppressMessages(library(optparse))
suppressMessages(library(ggrepel))
suppressMessages(library(harmony))

StopWatchEnd$LoadPackages  <- Sys.time()


#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--seuratObj",
                                type = "character",
                                default = NULL,
                                help = "path to seurat object to use as input (post-QC doublets)",
                                metavar= "character"
                               ),
                     make_option("--fileName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files, type NONE is you want to use sampleID",
                                metavar= "character"
                               ),
                     make_option("--outputDir",
                                type = "character",
                                help = "path to dir where all pipeline outputs will go",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--gene.filtering",
                                type = "double",
                                help = "percent to be removed fro mavg library size ie. 0.005",
                                default = NULL,
                                metavar= "double"
                               ),
                    make_option("--normalizationMethod",
                                type = "character",
                                help = "LogNormalize, scran or scTransform",
                                default = NULL,
                                metavar= "character"
                               ),
                     make_option("--varsRegress",
                                type = "character",
                                default = NULL,
                                help = "path to csv file with header REGRESS; or list NONE",
                                metavar= "character"
                               ),
                    make_option("--numVarGenes",
                                type = "integer",
                                default = NULL,
                                help = "number of varibale genes to identify",
                                metavar= "integer"
                               ),
                    make_option("--pcaScalingGenes",
                                type = "character",
                                help = "Use all or var for genes used in data scaling and PCA",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--computePC",
                                type = "integer",
                                help = "how many PCs to compute in PCA",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--numPC",
                                type = "integer",
                                help = "specify how many PCs to use for dim rediction, or list 0 for autodetection",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--addPC",
                                type = "integer",
                                help = "how many PCs to add to numPC to capture additional biological variation",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--outputPlots",
                                type = "logical",
                                help = "TRUE/FALSE to output plots",
                                default = NULL,
                                metavar= "logical"
                               ),
                    make_option("--minResolution",
                                type = "double",
                                help = "minimum clustering resolution; ie. 0.5",
                                default = NULL,
                                metavar= "double"
                               ),
                    make_option("--maxResolution",
                                type = "double",
                                default = NULL,
                                help = "maximum clustering resolution; ie. 1.0",
                                metavar= "double"
                               ),
                    make_option("--deGeneCutoff",
                                type = "integer",
                                default = NULL,
                                help = "cutoff for number DE genes per cluster; ie. 50",
                                metavar= "integer"
                               ),
                    make_option("--fdrCutoff",
                                type = "double",
                                default = NULL,
                                help = "Filter DE genes based on this FDR threshold, ir. 0.01",
                                metavar= "double"
                               ),
                   make_option("--optimalCluster",
                                type = "logical",
                                default = TRUE,
                                help = "TRUE/FALSE to select optimal clustering solution?",
                                metavar= "charater"
                                ),
                   make_option("--markerGenes",
                                type = "character",
                                default = NULL,
                                help = "Path to csv with marker genes, or use NONE",
                                metavar= "charater"
                               ),
                    make_option("--species",
                                type = "character",
                                default = NULL,
                                help = "species used in expeirment, current options include mouse or human)",
                                metavar= "charater"
                               )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

seuratObj <- opt$seuratObj
fileName <- opt$fileName
outputDir <- opt$outputDir
gene.filtering <- opt$gene.filtering
normalizationMethod <- opt$normalizationMethod
varsRegress <- opt$varsRegress
numVarGenes <- opt$numVarGenes
pcaScalingGenes <- opt$pcaScalingGenes
computePC <- opt$computePC
numPC  <- opt$numPC
addPC <- opt$addPC
outputPlots <- opt$outputPlots
minResolution <- opt$minResolution
maxResolution <- opt$maxResolution
deGeneCutoff <- opt$deGeneCutoff
fdrCutoff <- opt$fdrCutoff
optimalCluster <- opt$optimalCluster
markerGenes <- opt$markerGenes
species <- opt$species

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

print(paste0("Path to seurat object: ", seuratObj))

print(paste0("File prefix: ", fileName))

print(paste0("Output directory: ", outputDir))
print(paste0("Lowly Expressed Gene Cutoff: ", gene.filtering))
print(paste0("Normalization Method: ", normalizationMethod))


if (varsRegress == "NONE"){

    vars.regress <- NULL
    print(paste0("Variables to Regress: ", vars.regress))

} else {

    vars.regress <-  c(as.character(read.table(varsRegress, header = T)$REGRESS))
    print(paste0("Variables to Regress: ", paste(vars.regress, collapse = ", ")))

}


print(paste0("Number Variable Genes: ", numVarGenes))
print(paste0("Genes used for PCA and scaling: ", pcaScalingGenes))
print(paste0("Number of PCs to compute: ", computePC))
print(paste0("Number of PCs to use: ", ifelse(numPC == 0, "Auto-detect", numPC)))
print(paste0("Number of PCs to add: ", ifelse(normalizationMethod == "scTransform", "1.5x", addPC)))
print(paste0("Output Plots?: ", outputPlots))
print(paste0("Minimum Clustering Resolution: ", minResolution))
print(paste0("Maximum Clustering Resolution: ", maxResolution))
print(paste0("FDR Cutoff for Marker Genes: ", fdrCutoff))
print(paste0("Clusters must have a minimum of: ", deGeneCutoff, " DE genes"))
print(paste0("Run optimal clustering solution?: ", optimalCluster))

if (markerGenes == "NONE"){

    marker.genes <- NULL
    print(paste0("Assessing marker genes: ", markerGenes))

} else {

    marker.genes <- read.csv(markerGenes)
    print(paste0("Assessing marker genes: ", paste(marker.genes$GENE, collapse = ",")))

}

StopWatchEnd$SetUpVaribles  <- Sys.time()

setwd(outputDir)
dir.create(fileName)
setwd(fileName)
dir.create("data")
if(outputPlots == TRUE){dir.create("figures")}


#########################################
# Loads lineage-subsetted scRNA-seq data (seurat objects)
#########################################

print("")
print("*********************")
print("Read in Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadData <- Sys.time()


seurat.obj <- readRDS(seuratObj)

print("Total Dataset Size")
print(dim(seurat.obj@assays$RNA@data)) ## print out final size of object
  
StopWatchEnd$ReadData <- Sys.time()



#########################################
# Remove lowly expressed genes
#########################################

print("")
print("********************")
print("Filter out lowly expressed genes")
print(Sys.time())
print("********************")
print("")
StopWatchStart$FilterGenes <- Sys.time()

####remove gene.filtering % of lowest library size

nCells <- mean((table(seurat.obj@meta.data$Sample_ID))) * gene.filtering
print(paste0("Cutoff: ", gene.filtering, " of average cell size (", round(nCells, 2), " cells)"))

seurat.obj <- CreateSeuratObject(counts = seurat.obj@assays$RNA@counts,
                                 meta.data = seurat.obj@meta.data,
                                 min.cells = nCells
                                )

#seurat.obj@assays$RNA@counts <- seurat.obj@assays$RNA@counts[rowSums(as.matrix(seurat.obj@assays$RNA@counts > 0)) >= nCells, ]
#seurat.obj@assays$RNA@data <- seurat.obj@assays$RNA@data[rownames(seurat.obj@assays$RNA@counts), ]
print(dim(seurat.obj@assays$RNA@counts))
print(dim(seurat.obj@assays$RNA@data))

StopWatchEnd$FilterGenes <- Sys.time()


#########################################
# Normalization
#########################################

StopWatchStart$Normalization <- Sys.time()

if(normalizationMethod == "scran"){


    print("")
    print("********************")
    print("Normalize Data with Scran")
    print(Sys.time())
    print("********************")
    print("")

    # convert seurat counts to SingleCellExperiment object
    sce <- SingleCellExperiment(assays = list(counts = as.matrix(x = seurat.obj@assays$RNA@counts)))

    # quick cluster first; scran noramlizes within size factor pools
    # minClusterSize <- 0.05*ncol(seurat.obj@assays$RNA@counts)
    # print(paste0("Minimum Cluster Size for quickCluster: ", minClusterSize))
    minClusterSize <- 25
    print(paste0("Minimum Cluster Size for quickCluster: ", 25))

    clusters <- quickCluster(sce,
                             method = "igraph",
                             min.size= minClusterSize #hard coded to 25
                        )
    print(paste0("Identified ", length(unique(clusters)), " clusters..."))

    print("Computing Sum Factors...")
    sce <- computeSumFactors(sce, cluster=clusters)

    print("Normalizing...")
    #sce <- scater::normalize(sce, return_log = FALSE)
    # without log transform
    sce <- scater::logNormCounts(sce, log = FALSE)

    print("Update Master Seurat Object....")
    #add scran normalized counts into seurat obj
    seurat.obj <- NormalizeData(object = seurat.obj,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000
                               )

    #backup Seurat's and scran (pre-log) norm data
    seurat.obj@misc[["seurat_norm_data"]] <- as.matrix(x = seurat.obj@assays$RNA@data)
    seurat.obj@misc[["scran_norm_data_noLog"]] <- as.matrix(assay(sce, "normcounts"))

    #replace normalized data with scran
    seurat.obj@assays$RNA@data <- log(x = assay(sce, "normcounts") + 1)
    print(paste0("Consistency Check....",
                 identical(rownames(seurat.obj@assays$RNA@data),
                          rownames(seurat.obj@assays$RNA@counts))
         ))

    #add size factors from scran into metadata
    sizeFactors <- data.frame(sizeFactors(sce))
    colnames(sizeFactors) <- "scran_SizeFactor"
    rownames(sizeFactors) <- colnames(counts(sce))
    seurat.obj <- AddMetaData(seurat.obj, metadata = sizeFactors)

    #remove scran intermediate files to free up mem
    rm(sce)
    rm(clusters)
    rm(sizeFactors)

} else if(normalizationMethod == "LogNormalize"){


    print("")
    print("********************")
    print("Normalize Data with Seurat LogNormalize")
    print(Sys.time())
    print("********************")
    print("")

    seurat.obj <- NormalizeData(object = seurat.obj,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000
                               )


} else if(normalizationMethod == "scTransform"){


    print("")
    print("********************")
    print("Normalize Data with Seurat scTransform")
    print(Sys.time())
    print("********************")
    print("")
    print(paste0("Regressing out....", paste(vars.regress, collapse = ", ")))
    print(paste0("Using fixed number of variable genes... ", numVarGenes))
    seurat.obj <- SCTransform(seurat.obj,
                              new.assay.name = "SCT",
                              vars.to.regress = vars.regress,
                              variable.features.n = numVarGenes,
                              verbose = TRUE
                              #variable.features.rv.th = 1.3
                             )
}


StopWatchEnd$Normalization <- Sys.time()


#########################################
# VarGenes and Data Scaling
#########################################

if (normalizationMethod == "LogNormalize" | normalizationMethod == "scran"){

    print("")
    print("********************")
    print("Identify Variable Genes (vst)")
    print(Sys.time())
    print("********************")
    print("")
    StopWatchStart$VarGenes <- Sys.time()
    #use default Seurat Method, vst
    print(paste0("Calculating fixed number of variable genes... ", numVarGenes))
    seurat.obj <- FindVariableFeatures(seurat.obj,
                                       selection.method = "vst",
                                       nfeatures = numVarGenes
                                      )
    StopWatchEnd$VarGenes <- Sys.time()

    print("")
    print("********************")
    print("Scale Data")
    print(Sys.time())
    print("********************")
    print("")
    StopWatchStart$ScaleData <- Sys.time()
    print(paste0("Regressing out....", paste(vars.regress, collapse = ", ")))
    print(paste0("Genes included in scaling... ", pcaScalingGenes))


    if (pcaScalingGenes == "all"){

        all.genes <- rownames(seurat.obj)
        seurat.obj <- ScaleData(seurat.obj,
                                vars.to.regress = vars.regress,
                                features = all.genes
                               )

    } else if (pcaScalingGenes == "var"){

            seurat.obj <- ScaleData(seurat.obj,
                                vars.to.regress = vars.regress,
                                features =  VariableFeatures(object = seurat.obj)
                               )
    }


}

StopWatchEnd$ScaleData <- Sys.time()


#########################################
# Cell Cycle Scoring
#########################################

print("")
print("********************")
print("Score Cells for Cell Cycle Markers")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ScoreCellCycle <- Sys.time()
#As of now, there is no option to regress out cell cycle in this pipeline

# Also read in a list of cell cycle markers, from Tirosh et al, 2015

if (species == "human") {
 
  cc.genes <- readLines(con = "/cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/regev_lab_cell_cycle_genes.txt")
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]

}

if (species == "mouse") {
 
  cc.genes.human <- readLines(con = "/cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/regev_lab_cell_cycle_genes.txt")
  s.genes.human <- cc.genes.human[1:43]
  g2m.genes.human <- cc.genes.human[44:97]
  
  cc.genes <- read.table("/cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/regev_lab_cell_cycle_genes_both.txt")
  
  s.genes <- cc.genes[cc.genes$HGNC.symbol %in% s.genes.human, ]
  s.genes <- s.genes$MGI.symbol

  g2m.genes <- cc.genes[cc.genes$HGNC.symbol %in% g2m.genes.human, ]
  g2m.genes <- g2m.genes$MGI.symbol

  
}

seurat.obj <- CellCycleScoring(seurat.obj,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = FALSE
                              )
#calculate CC.Differennc
#just in case you want to regress with alt method in future
seurat.obj$CC.Difference <- seurat.obj$S.Score - seurat.obj$G2M.Score

StopWatchEnd$ScoreCellCycle <- Sys.time()




#########################################
# Run PCA
#########################################

print("")
print("********************")
print("Run PCA")
print(Sys.time())
print("********************")
print("")
StopWatchStart$PCA <- Sys.time()

print(paste0("Calculating... ", computePC, " PCs"))
#print(paste0("Calculating... ", 100, " PCs"))
print(paste0("Genes included in PCA... ", pcaScalingGenes))

if (pcaScalingGenes == "all"){

    all.genes <- rownames(seurat.obj)
    seurat.obj <- RunPCA(seurat.obj,
                         features = all.genes,
                         npcs = computePC,
                         #npcs = 100,
                         verbose = FALSE
                        )

} else if (pcaScalingGenes == "var"){

    seurat.obj <- RunPCA(seurat.obj,
                         features = VariableFeatures(object = seurat.obj),
                         npcs = computePC,
                         #npcs = 100,
                         verbose = FALSE
                        )
}


#########################################
# Determine significant PCs
#########################################


### isolate PCA data for output
pca <- data.frame(seurat.obj@reductions$pca@stdev)
pca$PC <- seq(1:nrow(pca))
pca <- pca[ ,c(2,1)]
colnames(pca)[2] <- "st.dev"
eigs <- pca$st.dev**2
pca$prop.var <- eigs / sum(eigs)


if(numPC == 0){

    print(paste0("Determining significant PCs...Auto (KneeArrower)"))


    cutoff.point <- findCutoff(pca$PC[1:computePC],
                               pca$st.dev[1:computePC],
                               method="first",
                               0.01 #derivative cutoff
                              )
    numPC <- round(cutoff.point$x)
    print(paste0(numPC, " significant PCs"))

} else if(numPC > 0){

    print(paste0("Determining significant PCs...User Defined"))
    print(paste0(numPC, " significant PCs"))
}


if(addPC > 0){

    print(paste0("Adding Additional PCs...", addPC, " PCs"))
    pc.use <- numPC + addPC
    print(paste0("PCs used for downstream analysis...", pc.use, " PCs"))

} else if (addPC == 0){

    pc.use <- numPC
    print(paste0("PCs used for downstream analysis...", pc.use, " PCs"))

}


if(normalizationMethod == "scTransform"){

    #add 1.5 times PCs for scTransform to capture more biology
    #Seurat recommendation
    pc.use <- round(numPC * 1.5)
    print(paste0("Multiply significant PCs 1.5x for scTransform..."))
    print(paste0("PCs used for downstream analysis...", pc.use, " PCs"))

}

StopWatchEnd$PCA <- Sys.time()


#########################################
# Batch correction (harmony)
#########################################

print("")
print("********************")
print("Running batch correction")
print(Sys.time())
print("********************")
print("")

# https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html
# http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/advanced.html

StopWatchStart$harmony <- Sys.time()
seurat.obj <- RunHarmony(seurat.obj, "SeqBatch", plot_convergence = F, theta = 0.5)
StopWatchEnd$harmony <- Sys.time()






#########################################
# tSNE and UMAP
#########################################

print("")
print("********************")
print("Non-linear Dimensionality Reduction")
print(Sys.time())
print("********************")
print("")

StopWatchStart$tSNE <- Sys.time()
print(paste0("Running tSNE..."))
seurat.obj <- RunTSNE(seurat.obj, reduction = "harmony", dims = 1:pc.use, verbose = FALSE, max_iter = 2000, check_duplicates = FALSE)
StopWatchEnd$tSNE <- Sys.time()

StopWatchStart$UMAP <- Sys.time()
print(paste0("Running UMAP..."))
seurat.obj <- RunUMAP(seurat.obj, reduction = "harmony", dims = 1:pc.use, verbose = FALSE)
StopWatchEnd$UMAP <- Sys.time()

print(paste0("Add coordinates to metadata..."))
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$umap@cell.embeddings)
                          )
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$tsne@cell.embeddings)
                          )
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$harmony@cell.embeddings)
                          )

#########################################
# Clustering
#########################################

print("")
print("********************")
print("Cluster Cells with Seurat")
print(Sys.time())
print("********************")
print("")
StopWatchStart$Clustering <- Sys.time()

# remove clustering results in meta.data from previous analyses
previous.clus.results <- colnames(seurat.obj@meta.data)[grep("_res", colnames(seurat.obj@meta.data))]
seurat.obj@meta.data <- seurat.obj@meta.data[, !colnames(seurat.obj@meta.data) %in% previous.clus.results]


#default in seurat is to use k = 20
print(paste0("Constructing SNN Graph with ", pc.use, " PCs..."))
seurat.obj <- FindNeighbors(seurat.obj,
                            reduction = "harmony",
                            dims = 1:pc.use,
                            verbose = TRUE
                           )

print(paste0("Find Clusters with ", pc.use, " PCs..."))
print(paste0("Clustering over a range of resolutions..."))


#Seurat recommends 0.4-1.2 resolution for 300k cells
#calculate 6 even ranges between to range specified by user

res.range <- round(seq(minResolution,
                       maxResolution,
                       length.out = 6
                      ),
                   2)
print(res.range)

seurat.obj <- FindClusters(seurat.obj,
                           resolution = res.range,
                           #method = "igraph", #better for large data
                           algorithm = 1, #default Louvian algorithm,
                           verbose = FALSE,
                           n.start = 50
                          )

##fix naming in meta.data

if (normalizationMethod == "LogNormalize" | normalizationMethod == "scran"){

    colnames(seurat.obj@meta.data) <- gsub("RNA_snn", "Seurat_cluster", colnames(seurat.obj@meta.data))
    seurat.obj@meta.data$seurat_clusters <- NULL

} else if (normalizationMethod == "scTransform"){

    colnames(seurat.obj@meta.data) <- gsub("SCT_snn", "Seurat_cluster", colnames(seurat.obj@meta.data))
    seurat.obj@meta.data$seurat_clusters <- NULL

}


StopWatchEnd$Clustering <- Sys.time()




#########################################
# DE Marker Genes
#########################################

if (optimalCluster == TRUE){
  
print("")
print("********************")
print("Differential Gene Expression Analysis")
print(Sys.time())
print("********************")
print("")
StopWatchStart$MarkerGenes <- Sys.time()

resolutions <- colnames(seurat.obj@meta.data)[grep("_res", colnames(seurat.obj@meta.data))]
#print(resolutions)

markers <- list()

for (i in 1:length(resolutions)){

    #set cell identity to match resolution
    seurat.obj <- SetIdent(seurat.obj, value = resolutions[i])
    print(paste0(i,
                 "/",
                 length(resolutions),
                 "...",
                 resolutions[i],
                 " (",
                 #length(unique(seurat.obj@meta.data[ ,resolutions[i]])),
                 length(unique(Idents(seurat.obj))),
                 " clusters)"
                ))

    DE.outs <- FindAllMarkers(seurat.obj,
                              test.use = "wilcox",
                              only.pos = TRUE,
                              min.pct = 0.25, #only test genes in 25% of both populations
                              logfc.threshold = 0.25, #log-fold change threshold
                              verbose = FALSE,
                              pseudocount.use = 1/nrow(seurat.obj@meta.data),
                              return.thresh = 0.01 #pvalue to return genes
                             )

    DE.outs$resolution <- as.numeric(gsub("Seurat_cluster_res.", "", resolutions[i]))
    DE.outs$sample <- fileName
    rownames(DE.outs) <- NULL
    markers[[resolutions[i]]] <- DE.outs
    rm(DE.outs) #help to conserve memory
    gc()

}

#combine gene marker lists from all solutions into one file
#save as RDS and csv
print("Merging DE results....")
markers <- do.call("rbind", markers)
rownames(markers) <- NULL

#remove genes that are above FDR cutoff
print(paste0("Filtering DE genes....FDR=", fdrCutoff))
markers_fdr <- markers[markers$p_val_adj <= fdrCutoff, ]
print(paste0("Removed ", nrow(markers) - nrow(markers_fdr), " genes...."))

StopWatchEnd$MarkerGenes <- Sys.time()



} else if (optimalCluster == FALSE){

  #if you are not picking optimal cluster solution
  #use the middle clustering solution within range as the one for plotting
  middle.res <- res.range[length(res.range)/2]
  print(paste0("DE testing not required if optimal clustering not performed"))

}

###### HERE IS WHERE WE NEED TO DETERMIN IF LOWER RES NEEDS
###### TO BE CALCULATED
###### TABL

#########################################
# Sil Width
#########################################

if (optimalCluster == TRUE){

print("")
print("********************")
print("Calculate Silhouette Width")
print("********************")
print(Sys.time())
print("")
StopWatchStart$SilhouetteWidth <- Sys.time()
#?silhouette() in cluster package
#this will give a value for each cell
#https://rpkgs.datanovia.com/factoextra/reference/fviz_silhouette.html
#calculate sil values using pca matrix
#same matrix used for input to clustering

pca.dist <- as.matrix(seurat.obj@reductions$pca@cell.embeddings[ ,1:pc.use])
pca.dist <- dist(pca.dist)

silhouette.width <- list()

for(i in 1:length(resolutions)){

    print(paste0(i,
                 "/",
                 length(resolutions),
                 "...",
                 resolutions[i]
                ))

    cl <- seurat.obj@meta.data[ ,resolutions[i]] #beware this adds +1 to cluster numbers

    sil <- cluster::silhouette(as.integer(cl), pca.dist)
    sil <- data.frame(sil[, 1:ncol(sil)])
    sil$cluster <- sil$cluster - 1
    sil$resolution <- as.numeric(gsub("Seurat_cluster_res.", "", resolutions[i]))
    sil$cell <- rownames(pca.dist)

    silhouette.width[[resolutions[i]]] <- sil
    rm(sil) #help to conserve memory
    gc()

}

print("Merging silhouette widths across solutions....")
silhouette.width <- do.call("rbind", silhouette.width)
rownames(silhouette.width) <- NULL
head(silhouette.width)

StopWatchEnd$SilhouetteWidth <- Sys.time()



#########################################
# Select Clustering Solution
#########################################

print("")
print("********************")
print("Select Optimal Clustering Solution")
print(Sys.time())
print("********************")
print("")
StopWatchStart$OptimalSolution <- Sys.time()

print(paste0("Minimum ",
             deGeneCutoff,
             " DE genes/cluster after FDR filtering....")
     )

#extract the number of DE genes per cluster
deGenes <- table(markers_fdr$cluster, markers_fdr$resolution)
deGenes <-  reshape2::melt(deGenes)
colnames(deGenes) <- c("cluster", "resolution", "nDEGenes")

#filter on cluster number

ress <- unique(markers_fdr$resolution)
newdata <- list()
for (i in 1:length(ress)){

    clust.name <- paste0("Seurat_cluster_res.", ress[i])
    #print(clust.name)
    #subset <- seurat.obj@meta.data[ ,grep(clust.name, colnames(seurat.obj@meta.data))]
    subset <- seurat.obj@meta.data[ ,clust.name]
    num.clust <- length(unique(subset))
    #print(num.clust)
    a <- deGenes[deGenes$resolution == ress[i], ]
    a <- a[1:num.clust, ]
    newdata[[i]] <- a

}

deGenes <- do.call(rbind, newdata)

#deGenes <- deGenes %>% drop_na()
#colnames(deGenes) <- c("cluster", "resolution", "nDEGenes")

print("Cluster for each resolution with lowest # of DE genes...")
clust.min <- data.frame(deGenes %>% group_by(resolution) %>% top_n(n = -1, wt = nDEGenes))
clust.min$Pass_DE <- clust.min$nDEGenes >= deGeneCutoff
clust.min <- clust.min[!duplicated(clust.min$resolution), ]
print(dim(clust.min))
print("")

possible.solutions <- clust.min$resolution[clust.min$nDEGenes >= deGeneCutoff]
print(paste0("Resolutions that pass # DE gene cutoff....", paste(possible.solutions, collapse = ", ")))

#print(paste0("Calculate min and max....", paste(possible.solutions, collapse = ", ")))
print("Calculate median silhuoette width for resolutions....")

silhouette.width$res_cluster <- paste0("res",
                                       silhouette.width$resolution,
                                       "_",
                                       silhouette.width$cluster
                                      )

##first calculate sil width for each cluster
#sil_aggr <- aggregate(silhouette.width[, "sil_width"],
#                      list(silhouette.width$res_cluster),
#                      mean
#                     )
#colnames(sil_aggr) <- c("res_cluster", "sil_avg")
#sil_aggr$resolution <- sapply(strsplit(sil_aggr$res_cluster,"_"), `[`, 1)
#sil_aggr <- sil_aggr[order(sil_aggr$resolution), ]

# then average the sil width across all clusters
# this gives cluster even weight in the calculation
#sil_avg <- aggregate(sil_aggr[, "sil_avg"],
#                      list(sil_aggr$resolution),
#                      mean
#                     )

#clust.min$Sil_Avg <- sil_avg$x

#### MEDIAN CALCULATIONS

sil_aggr <- aggregate(silhouette.width[, "sil_width"],
                      list(silhouette.width$res_cluster),
                      median
                     )
colnames(sil_aggr) <- c("res_cluster", "sil_median")
sil_aggr$resolution <- sapply(strsplit(sil_aggr$res_cluster,"_"), `[`, 1)
sil_aggr <- sil_aggr[order(sil_aggr$resolution), ]

sil_med <- aggregate(sil_aggr[, "sil_median"],
                      list(sil_aggr$resolution),
                      median
                     )
clust.min$Sil_median <- sil_med$x

#### PICK OPTIMAL SOLUTION
print(clust.min)
opt.solution <- clust.min[clust.min$Pass_DE == TRUE, ]
opt.solution <- opt.solution[opt.solution$Sil_median == max(opt.solution$Sil_median),  ] #has to pass DE test

# if multiple solutions are tied...
# pick the highest resolution
if(nrow(opt.solution) > 1){
    opt.solution <- opt.solution[opt.solution$resolution == max(opt.solution$resolution), ]
}

#add optimal resolution to metadata and assign ident
seurat.obj@meta.data$Optimal_res <- opt.solution$resolution
new.ident <- colnames(seurat.obj@meta.data)[grep(opt.solution$resolution, colnames(seurat.obj@meta.data))]
seurat.obj <- SetIdent(seurat.obj, value = new.ident)

print(paste0("Optimal solution is.....Resolution = ", opt.solution$resolution))

deGenes$res_cluster <- paste0("res", deGenes$resolution, "_", deGenes$cluster)
deGenes <- deGenes[order(deGenes$res_cluster), ]
deGenes$sil_median <- sil_aggr$sil_median
optimal.file <- paste0("./data/",fileName, "_DE_Silhouette.csv")
print(optimal.file)
write.csv(deGenes, file = optimal.file)

StopWatchEnd$OptimalSolution <- Sys.time()

} else if (optimalCluster == FALSE){

  #if you are not picking optimal cluster solution
  #use the middle clustering solution within range as the one for plotting
  middle.res <- res.range[length(res.range)/2]
  print(paste0("Picking middle of resolution range as optimal clustering solution: ", middle.res))
  seurat.obj@meta.data$Optimal_res <- middle.res

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
raw <- Matrix(seurat.obj@assays$RNA@counts, sparse = TRUE)
write10xCounts(path = raw.count.dir,
               x = raw
               )

print("Normalized Counts....")
norm.count.dir <- paste0("./data/", fileName, "_normCounts")
norm.count <- Matrix(seurat.obj@assays$RNA@data, sparse = TRUE)
write10xCounts(path = norm.count.dir,
               x = norm.count
               )

print("Meta Data....")
meta <- data.frame(seurat.obj@meta.data)
meta.file <- paste0("./data/", fileName, "_metaData.csv")
write.csv(meta, file = meta.file)
meta.file.2 <- paste0("./data/", fileName, "_metaData.rds")
saveRDS(meta, file = meta.file.2)

print("Scree Plot Information....")
scree.file <- paste0("./data/", fileName, "_scree.csv")
write.csv(pca, file = scree.file)

print("PCA....")
pca.file <- paste0("./data/", fileName, "_PCA.rds")
pca.mat <- seurat.obj@reductions$pca
saveRDS(pca.mat, file = pca.file)

if (optimalCluster == TRUE){
print("DE Markers....")
DE.file <- paste0("./data/", fileName, "_DE_markers.csv")
write.csv(markers, file = DE.file)
DE.file.2 <- paste0("./data/", fileName, "_DE_markers_FDR_filered.csv")
write.csv(markers_fdr, file = DE.file.2)


print("Silhouette Width....")
sil.file <- paste0("./data/", fileName, "_silhouette_widths.rds")
saveRDS(silhouette.width, file = sil.file)
sil.file.2 <- paste0("./data/", fileName, "_clustering_dashboard.csv")
write.csv(deGenes, file = sil.file.2)
}

print("Saving Seurat Object....")
seurat.file <- paste0("./data/", fileName, "_seurat.rds")
saveRDS(seurat.obj, file = seurat.file)

StopWatchEnd$SaveData <- Sys.time()



#########################################
# Plotting
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


    if(class(marker.genes) == "data.frame"){

        print("Adding Marker Genes to Metadata...")
        gene.index <- as.character(marker.genes$GENE) %in% rownames(seurat.obj@assays$RNA@data)
        genes <- as.character(marker.genes$GENE)[gene.index]
        subset <- seurat.obj@assays$RNA@data[genes, ]
        subset <- data.frame(data.matrix(t(subset)))
        seurat.obj <- AddMetaData(seurat.obj,
                                  metadata = subset
                                 )

        if(FALSE %in% gene.index){
            print("The following genes are not in the dataset: ")
            print(as.character(marker.genes$GENE)[!gene.index])
        }



    }

    meta <- data.frame(seurat.obj@meta.data)

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

    #######################

    print("2/6....Variable Gene Plot")

    top20 <- head(VariableFeatures(seurat.obj), 20)

    VAR1 <- VariableFeaturePlot(seurat.obj)
    VAR2 <- LabelPoints(plot = VAR1, points = top20, repel = TRUE)

    var.pdf <- paste0("./figures/", fileName, "_VarGenes.pdf")
    print(var.pdf)
    pdf(var.pdf, width = 10, height = 8)
    print(VAR2)
    dev.off()

    #######################

    print("3/6....PCA Scree Plot")

    pc.pdf <- paste0("./figures/", fileName, "_PCA_plots.pdf")
    print(pc.pdf)
    pdf(pc.pdf, width = 10, height = 5)
    par(mfrow=c(1,2))

            plot(pca$PC,
                 pca$st.dev,
                 xlab = "Principal Components",
                 ylab = "Standard Deviation",
                 main = fileName
                    )
            abline(v=round(pc.use), lty = 2)
            legend("topright",
                   legend = c(paste0(round(pc.use), " PCs")),
                   bty = 'n')

        plot(pca$PC,
             pca$prop.var,
             xlab = "Principal Components",
             ylab = "Proportion Variance"
            )
        abline(v=round(pc.use), lty = 2)
        legend("topright",
               legend = c(paste0(round(sum(pca$prop.var[1:round(pc.use)]), 3), " cumulative variance")),
               bty = 'n')

    dev.off()



    #######################
    if (optimalCluster == TRUE){

    print("4/6....Optimal Clustering Dashboard")

    nClusters <- c()
    resolutions <- unique(deGenes$resolution)

    cols <- ifelse(resolutions == unique(meta$Optimal_res), "#7fcdbb", "grey")

    for (i in 1:length(unique(deGenes$resolution))){
        nClusters[i] <- nrow(deGenes[deGenes$resolution == resolutions[i], ])
    }

    plot.dat <- data.frame(nClusters, resolutions)
    plot.dat$resolutions <- factor(plot.dat$resolutions)

    dashboard.1 <- ggplot(plot.dat, aes (x=resolutions, y = nClusters)) +
                    geom_bar(stat = "identity", color = "black", fill = cols) + theme_classic() +
                     geom_text(aes(label=nClusters), vjust=1.6, color = "black", size=3.5)+
                    xlab("") + ylab("Clusters") +
                    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
                    ggtitle(fileName)

    ##### boxplot of DE genes
    deGenes$resolution <- factor(deGenes$resolution)

    dashboard.2 <- ggplot(deGenes, aes (x=resolution, y = nDEGenes)) +
                   geom_boxplot(outlier.shape = NA, fill = cols) + theme_classic() +
                 #geom_text(aes(label=nClusters), vjust=1.6, color = "white", size=3.5)+
                    xlab("") + ylab("DE Genes / Cluster") +
                    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
                    geom_jitter(shape=16, position=position_jitter(0.2)) +
                    geom_hline(yintercept = deGeneCutoff, col = "red", lty = 2)


    dashboard.3 <- ggplot(deGenes, aes (x=resolution, y = sil_median)) +
                   geom_boxplot(outlier.shape = NA, fill = cols) + theme_classic() +
                     #geom_text(aes(label=nClusters), vjust=1.6, color = "white", size=3.5)+
                    xlab("Resolutions") +
                    ylab("Median Silhouette Width / Cluster") +
                    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
                    geom_jitter(shape=16, position=position_jitter(0.2))
                    #geom_hline(yintercept = deGeneCutoff, col = "red", lty = 2)

    dashboard.pdf <- paste0("./figures/", fileName, "_ClusterResolution_Dashboard.pdf")
    print(dashboard.pdf)

    pdf(dashboard.pdf, width = 10, height = 10)

    print(ggarrange(dashboard.1,
              dashboard.2,
              dashboard.3,
              ncol = 1,
              nrow = 3
             ))

    dev.off()

    } else if (optimalCluster == FALSE){

      print("4/6....Skip Plotting Optimal Clustering Dashboard")

    }
    #######################

    print("5/6.....tSNEs and UMAPs")


  #### PLOT tSNEs

    meta <- seurat.obj@meta.data
    meta$cluster <- meta[ ,colnames(meta) == paste0("Seurat_cluster_res.", unique(meta$Optimal_res))]
    clust.col <- colnames(meta)[grep(unique(meta$cluster), colnames(meta))]
    cent <- meta %>% group_by(cluster) %>% select(tSNE_1, tSNE_2) %>% summarize_all(median)
    
  tsne.c <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "cluster")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Clusters") + ggtitle(paste0(fileName, " - Cluster: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) +
            geom_label_repel(aes(label = cluster), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none') 
  
  tsne.s <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Sample_ID")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Samples") + ggtitle(paste0(fileName, " - Sample ID: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
            #geom_label_repel(aes(label = Sample_ID), data = cent, show.legend = F, size = 2)  +
            #theme(legend.position='none') 
  
  tsne.d <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Cohort_ID")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Disease Group") + ggtitle(paste0(fileName, " - Disease Group: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
            #geom_label_repel(aes(label = Cohort_ID), data = cent, show.legend = F, size = 2)  +
            #theme(legend.position='none') 

  tsne.pdf <- paste0("./figures/", fileName, "_tSNE.pdf")
  pdf(tsne.pdf, width = 10, height = 10)
  print(tsne.c)
  print(tsne.s)
  print(tsne.d)
  dev.off()

  
  
  
  #### PLOT UMAPs
    
    cent <- meta %>% group_by(cluster) %>% select(UMAP_1, UMAP_2) %>% summarize_all(median)
    
  umap.c <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "cluster")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Clusters") + ggtitle(paste0(fileName, " - Cluster: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) +
            geom_label_repel(aes(label = cluster), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none') 
  
  umap.s <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Sample_ID")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Samples") + ggtitle(paste0(fileName, " - Sample ID: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
            #geom_label_repel(aes(label = Sample_ID), data = cent, show.legend = F, size = 2)  +
            #theme(legend.position='none') 
  
 umap.d <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Cohort_ID")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Disease Group") + ggtitle(paste0(fileName, " - Disease Group: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
            #geom_label_repel(aes(label = Cohort_ID), data = cent, show.legend = F, size = 2)  +
            #theme(legend.position='none') 

  umap.pdf <- paste0("./figures/", fileName, "_UMAP.pdf")
  pdf(umap.pdf, width = 10, height = 10)
  print(umap.c)
  print(umap.s)
  print(umap.d)  
  dev.off()
  
  
  
    tsne.tiff <- paste0("./figures/", fileName, "_Cluster_tSNE.tiff")
    tiff(tsne.tiff, width = 700, height = 600)
    print(tsne.c)
    dev.off()
  
    tsne.tiff <- paste0("./figures/", fileName, "_Sample_tSNE.tiff")
    tiff(tsne.tiff, width = 700, height = 600)
    print(tsne.s)
    dev.off()

    tsne.tiff <- paste0("./figures/", fileName, "_Cohort_tSNE.tiff")
    tiff(tsne.tiff, width = 700, height = 600)
    print(tsne.d)
    dev.off()

  
    umap.tiff <- paste0("./figures/", fileName, "_Cluster_umap.tiff")
    tiff(umap.tiff, width = 700, height = 600)
    print(umap.c)
    dev.off()
  
    umap.tiff <- paste0("./figures/", fileName, "_Sample_umap.tiff")
    tiff(umap.tiff, width = 700, height = 600)
    print(umap.s)
    dev.off()

    umap.tiff <- paste0("./figures/", fileName, "_Cohort_umap.tiff")
    tiff(umap.tiff, width = 700, height = 600)
    print(umap.d)
    dev.off()
  
  
  #######################


    if(class(marker.genes) == "data.frame"){

        #######################

        print("6/6.....Marker Gene Expression")

        sigs <- colnames(meta)[colnames(meta) %in% marker.genes$GENE]
        print(sigs)

        #### UMAP with marker gene expression in hexbin
        markers.pdf <- paste0("./figures/", fileName, "_GeneMarkers_UMAP.pdf")
        print(markers.pdf)
        pdf(markers.pdf, height = 5, width = 5)

        for (i in 1:length(sigs)){

                plot.title <-  paste(sigs[i], "\n", fileName, "  ", nrow(meta), "cells")
                p <-  ggplot(meta, aes_string(x="UMAP_1", y="UMAP_2", z=sigs[i])) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = sigs[i])
                print(p)

            }

        dev.off()

       #### tSNE with marker gene expression in hexbin
      markers.pdf.t <- paste0("./figures/", fileName, "_GeneMarkers_tSNE.pdf")
      print(markers.pdf.t)
      pdf(markers.pdf.t, height = 5, width = 5)

        for (i in 1:length(sigs)){

                plot.title <-  paste(sigs[i], "\n", fileName, "  ", nrow(meta), "cells")
                p <-  ggplot(meta, aes_string(x="tSNE_1", y="tSNE_2", z=sigs[i])) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = sigs[i])
                print(p)

            }

    dev.off()

    StopWatchEnd$Plotting <- Sys.time()

    }

}


#### Tiff files of tSNE with marker gene expression in hexbin

setwd("figures")
dir.create("markers")

for (i in seq_along(sigs)) {
    
    goi <- sigs[i] 
    
    meta_ordered <- meta[order(meta[,goi], decreasing = F), ]

    plot.title <-  paste(goi, "\n", fileName, "  ", nrow(meta), "cells")
    
    p <-  ggplot(meta_ordered, aes_string(x="tSNE_1", y="tSNE_2", z=goi)) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = goi)
                print(p)
    
    mypath <- paste0("./markers/", fileName, "_", goi, "_GeneMarkers_tSNE.tiff")

    tiff(file=mypath)
    print(p)
    dev.off()
    
}


#### Tiff files of umap with marker gene expression in hexbin

for (i in seq_along(sigs)) {
    
    goi <- sigs[i] 
    
    meta_ordered <- meta[order(meta[,goi], decreasing = F), ]

    plot.title <-  paste(goi, "\n", fileName, "  ", nrow(meta), "cells")
    
    p <-  ggplot(meta_ordered, aes_string(x="UMAP_1", y="UMAP_2", z=goi)) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = goi)
                print(p)
    
    mypath <- paste0("./markers/", fileName, "_", goi, "_GeneMarkers_UMAP.tiff")

    tiff(file=mypath)
    print(p)
    dev.off()
    
  }




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

