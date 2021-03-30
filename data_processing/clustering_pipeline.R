###########################################################
#             Script to run clustering pipeline           #
#                L.Richards (D.Croucher mods)             #
#                       Sept 2020                         #
###########################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads pre-QC'ed scRNA-seq data (seurat objects) and filters lowly expressed genes
### 2) Normalizes (LogNormalize, scran, scTransform) and scales expression matrix
### 3) Identifies variable genes (vst)
### 4) Cell cyle scoring
### 5) Run PCA and determine number of significant PCs (can auto-detect scree)
### 6) Non-linear dimensionality reduction (tSNE, UMAP)
### 6) Cluster cells over range of 10 resolutions
### 7) Runs differential gene expression (DGE) over all clustering resolutions
### 8) Calculates silhouette width for clusters
### 9) Selects optimal clustering resolution based on DE-gene cutoff and max silhouette width
### 10) Output plots
### 11) Saves data
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
## Rscript /cluster/projects/pughlab/projects/scVKMYC/scripts/clustering_pipeline.R --samples /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/vkmyc_cohort_sampleInput.txt \
##     --fileName scVKMYC_Cohort \
##     --inputPath /cluster/projects/pughlab/projects/scVKMYC/analysis/test/seurat/scVKMYC_test_preprocess/data/ \
##     --gene.filtering 0.001 \
##     --normalizationMethod LogNormalize \
##     --varsRegress /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/varsRegress.txt \
##     --numVarGenes 3000 \
##     --pcaScalingGenes var \
##     --computePC 75 \
##     --numPC 0 \
##     --addPC 0 \
##     --outputPlots FALSE \
##     --minResolution 1 \
##     --maxResolution 3.5 \
##     --deGeneCutoff 15 \
##     --fdrCutoff 0.05 \
##     --outputDir /cluster/projects/pughlab/projects/scVKMYC/analysis/test/seurat/ \
##     --markerGenes /cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/CellTypeMarkerGenes_AllCells.csv \
##     --optimalCluster FALSE \
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

StopWatchEnd$LoadPackages  <- Sys.time()


#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--samples",
                                type = "character",
                                default = NULL,
                                help = "path to csv with column SAMPLEID; each sample on own line",
                                metavar= "character"
                               ),
                     make_option("--fileName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files, type NONE is you want to use sampleID",
                                metavar= "character"
                               ),
                     make_option("--inputPath",
                                type = "character",
                                help = "path to dir with seurat objects to read in",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--gene.filtering",
                                type = "double",
                                help = "percent to be removed from avg library size ie. 0.005",
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
                    make_option("--outputDir",
                                type = "character",
                                help = "path to dir where all pipeline outputs will go",
                                default = NULL,
                                metavar= "character"
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
                               ),
                     make_option("--optimalCluster",
                                type = "logical",
                                default = TRUE,
                                help = "TRUE/FALSE to select optimal clustering solution?",
                                metavar= "charater"
                                )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

samples <- opt$samples
fileName <- opt$fileName
inputPath <- opt$inputPath
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
markerGenes <- opt$markerGenes
optimalCluster <- opt$optimalCluster
outputDir <- opt$outputDir
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

sampleID <- c(as.character(read.table(samples, header = T)$SAMPLEID))
if(length(sampleID) == 1){
    print(paste0("Sample IDs: ", sampleID))
} else if(length(sampleID > 1)){
    print(paste0("Sample IDs: ", paste(sampleID, collapse = ", ")))
}

if(fileName == "NONE"){assign("fileName", sampleID)}
print(paste0("File prefix: ", fileName))

print(paste0("Input Path: ", inputPath))
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
#print(paste0("Run optimal clustering solution?: ", optimalCluster))

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
# Read in Data
#########################################

print("")
print("*********************")
print("Read in Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadData <- Sys.time()

files <- list.files(inputPath, pattern = "seurat.rds")
seurat.dir <- paste0(inputPath, files)
seurat.obj <- readRDS(seurat.dir)

print("Total Dataset Size")
print(dim(seurat.obj@assays$RNA@data)) ## print out final size of object

StopWatchEnd$ReadData <- Sys.time()


#########################################
# Remove QC Failed Cells
#########################################

print("")
print("********************")
print("Remove Failed QC Cells")
print(Sys.time())
print("********************")
print("")
StopWatchStart$RemoveFlaggedCells <- Sys.time()

print(table(seurat.obj@meta.data$Cell_QC))
seurat.obj <- subset(seurat.obj , subset = Cell_QC == "PASS")
print(dim(seurat.obj@assays$RNA@data))

StopWatchEnd$RemoveFlaggedCells <- Sys.time()


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

#print(paste0("Calculating... ", computePC, " PCs"))
print(paste0("Calculating... ", 200, " PCs"))
print(paste0("Genes included in PCA... ", pcaScalingGenes))

if (pcaScalingGenes == "all"){

    all.genes <- rownames(seurat.obj)
    seurat.obj <- RunPCA(seurat.obj,
                         features = all.genes,
                         #npcs = computePC,
                         npcs = 200,
                         verbose = FALSE
                        )

} else if (pcaScalingGenes == "var"){

    seurat.obj <- RunPCA(seurat.obj,
                         features = VariableFeatures(object = seurat.obj),
                         #npcs = computePC,
                         npcs = 200,
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
seurat.obj <- RunTSNE(seurat.obj, dims = 1:pc.use, verbose = FALSE, check_duplicates = FALSE)
StopWatchEnd$tSNE <- Sys.time()

StopWatchStart$UMAP <- Sys.time()
print(paste0("Running UMAP..."))
seurat.obj <- RunUMAP(seurat.obj, dims = 1:pc.use, verbose = FALSE)
StopWatchEnd$UMAP <- Sys.time()

print(paste0("Add coordinates to metadata..."))
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$umap@cell.embeddings)
                          )
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$tsne@cell.embeddings)
                          )



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

    qc.plots <- paste0("./figures/",fileName, "_QC.tiff")
    tiff(qc.plots, width=720,height=400)
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


    #### Sample ID

    tSNE.1 <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Sample_ID")) +
          geom_point(alpha = 0.5) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs(color = "Sample_ID") +
          ggtitle(paste0(fileName, "TSNE Coloured by Sample ID", sep = " "))

    umap.1 <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Sample_ID")) +
          geom_point(alpha = 0.5) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs(color = "Sample_ID") +
          ggtitle(paste0(fileName, "UMAP Coloured by Sample ID", sep = " "))


    ## Cohort ID

    tSNE.2 <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Cohort_ID")) +
          geom_point(alpha = 0.5) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs(color = "Cohort_ID") +
          ggtitle(paste0(fileName, "TSNE Coloured by Disease State", sep = " "))

    umap.2 <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Cohort_ID")) +
          geom_point(alpha = 0.5) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs(color = "Cohort_ID") +
          ggtitle(paste0(fileName, "UMAP Coloured by Disease State", sep = " "))


    tsne.pdf <- paste0("./figures/", fileName, "_tSNE.pdf")
    pdf(tsne.pdf, width = 10, height = 10)
    print(tSNE.1)
    print(tSNE.2)
    dev.off()

    umap.pdf <- paste0("./figures/", fileName, "_UMAP.pdf")
    pdf(umap.pdf, width = 10, height = 10)
    print(umap.1)
    print(umap.2)
    dev.off()

  
  
  
    tsne.tiff <- paste0("./figures/", fileName, "_Sample_tSNE.tiff")
    tiff(tsne.tiff, width = 700, height = 600)
    print(tSNE.1)
    dev.off()

    tsne.tiff <- paste0("./figures/", fileName, "_Cohort_tSNE.tiff")
    tiff(tsne.tiff, width = 700, height = 600)
    print(tSNE.2)
    dev.off()

  
    umap.tiff <- paste0("./figures/", fileName, "_Sample_umap.tiff")
    tiff(umap.tiff, width = 700, height = 600)
    print(umap.1)
    dev.off()

    umap.tiff <- paste0("./figures/", fileName, "_Cohort_umap.tiff")
    tiff(umap.tiff, width = 700, height = 600)
    print(umap.2)
    dev.off()
  
 
  
    #######################


    if(class(marker.genes) == "data.frame"){

    setwd("figures")
    dir.create("markers")
      
        #######################

        print("6/6.....Marker Gene Expression")

        sigs <- colnames(meta)[colnames(meta) %in% marker.genes$GENE]
        print(sigs)

        #### UMAP with marker gene expression in hexbin
        markers.pdf <- paste0("./markers/", fileName, "_GeneMarkers_UMAP.pdf")
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
      markers.pdf.t <- paste0("./markers/", fileName, "_GeneMarkers_tSNE.pdf")
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

         #### Tiff files of tSNE with marker gene expression in hexbin

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

setwd("../")

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

#print("DE Markers....")
#DE.file <- paste0("./data/", fileName, "_DE_markers.csv")
#write.csv(markers, file = DE.file)
#DE.file.2 <- paste0("./data/", fileName, "_DE_markers_FDR_filered.csv")
#write.csv(markers_fdr, file = DE.file.2)

if (optimalCluster == TRUE){
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
