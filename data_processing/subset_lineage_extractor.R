###########################################################
#  Script to extract cell lineages for subset analysis    #
#           D.Croucher (L.Richards template)              #
#                      Sept 2020                          #
###########################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load singlet QC'ed and SingleR annotated scRNA-seq data and integrate lineage subset - cluster assignments 
### 2) Plot any new marker genes
### 3) Subset cohort object into lineage-specific objects 
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
## Rscript /cluster/projects/pughlab/projects/scVKMYC/scripts/subset_lineage_extractor.R \
## --seuratObj /cluster/projects/pughlab/projects/scVKMYC/analysis/test/singler/singlecell_mode/scVKMYC_test/scVKMYC_test.ImmGenData_SinglerAnnotations_seurat.rds \
## --lineageSubsets /cluster/projects/pughlab/projects/scVKMYC/analysis/test/vkmyc_test_lineageSubsets.csv \
## --outName scVKMYC_test_lineages \
## --outputDir /cluster/projects/pughlab/projects/scVKMYC/analysis/test/lineages \
## --markerGenes /cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/CellTypeMarkerGenes_AllCells.csv \
## --updateMarkerGenes TRUE \
## --markerGeneDir /cluster/projects/pughlab/projects/scVKMYC/analysis/test/seurat/scVKMYC_test_multires/figures/
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

option_list <- list(make_option("--seuratObj",
                                type = "character",
                                default = NULL,
                                help = "path to final seurat object containing all cell types to use as input (post-QC singlets with multires clustering)",
                                metavar= "character"
                               ),
                    make_option("--lineageSubsets",
                                type = "character",
                                default = NULL,
                                help = "path to csv file containing cluster to cell type assignments",
                                metavar= "character"
                               ),
                    make_option("--outName",
                                type = "character",
                                help = "will be appended to all output files",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--outputDir",
                                type = "double",
                                help = "path to dir where all pipeline outputs will go",
                                default = NULL,
                                metavar= "double"
                               ),
                    make_option("--markerGenes",
                                type = "character",
                                default = NULL,
                                help = "Path to csv with marker genes, or use NONE",
                                metavar= "charater"
                               ),
                    make_option("--updateMarkerGenes",
                                type = "logical",
                                default = TRUE,
                                help = "TRUE/FALSE to output additional set of marker genes",
                                metavar= "charater"
                               ),
                    make_option("--markerGeneDir",
                                type = "character",
                                default = TRUE,
                                help = "path to dir where all marker gene feature plots will go",
                                metavar= "charater"
                               )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

seuratObj <- opt$seuratObj
lineageSubsets <- opt$lineageSubsets
outName <- opt$outName
outputDir <- opt$outputDir
markerGenes <- opt$markerGenes
updateMarkerGenes <- opt$updateMarkerGenes
markerGeneDir <- opt$markerGeneDir

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

print(paste0("Reading in lineage subsets from: ", lineageSubsets))

lineage.subset <- read.csv(lineageSubsets)
unique.lineages <- unique(lineage.subset$Lineage)

print(paste0("The following lineages will be extracted: ", paste(unique.lineages, collapse = ",")))

print(paste0("File prefix: ", outName))

print(paste0("Output directory: ", outputDir))

if (updateMarkerGenes == TRUE){
 
    marker.genes <- read.csv(markerGenes)
    markerGeneDir <- markerGeneDir
    print(paste0("The following markers genes will be saved to ", markerGeneDir))
    print("")
    print(paste(marker.genes$GENE, collapse = ","))
  
} else {

    print(paste0("No marker genes to update"))

}


StopWatchEnd$SetUpVaribles  <- Sys.time()

setwd(outputDir)
dir.create(outName)
setwd(outName)


#########################################
# Load singlet QC'ed and SingleR annotated scRNA-seq data and integrate lineage subset - cluster assignments 
#########################################

print("")
print("*********************")
print("Read in Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadData <- Sys.time()

seurat <- readRDS(seuratObj)

print(paste0("Total cells: ", ncol(seurat)))
print(paste0("Total genes: ", nrow(seurat)))


#lift lineage subset - cluster assignments over to each cell (row in meta.data)

meta <- seurat@meta.data
meta$cluster <- meta[ ,colnames(meta) == paste0("Seurat_cluster_res.", unique(meta$Optimal_res))]
meta$cluster <- as.numeric(as.character(meta$cluster)) + 1

clusters <- sort(unique(meta$cluster)) #vector of clusters ordered 1 to n+1

for (i in seq_along(clusters)) {
  
    meta[meta$cluster %in% clusters[i], "Subset_lineage"] <- lineage.subset[i,"Lineage"]

}

seurat <- AddMetaData(seurat, meta)
seurat@meta.data$cluster <- meta[ ,colnames(meta) == paste0("Seurat_cluster_res.", unique(meta$Optimal_res))]

print(table(seurat@meta.data$Subset_lineage))
print(table(seurat@meta.data$Subset_lineage, seurat@meta.data$Sample_ID))
print(prop.table(x = table(seurat@meta.data$Subset_lineage, seurat@meta.data$Sample_ID), margin = 2))

StopWatchEnd$ReadData <- Sys.time()


#########################################
# Plot any new marker genes
#########################################

if (updateMarkerGenes == FALSE) {

  print("")
  print("*********************")
  print("No new marker genes to plot")
  print("********************")
  print("") 
  
}else{

  print("")
  print("*********************")
  print("New marker genes indicated")
  print(Sys.time())
  print("********************")
  print("")
  StopWatchStart$PlotNewMarkers <- Sys.time()

########## Add marker gene expression values to meta data

  print("Adding Marker Genes to Metadata...")
  
  gene.index <- as.character(marker.genes$GENE) %in% rownames(seurat@assays$RNA@data)
  genes <- as.character(marker.genes$GENE)[gene.index]
  subset <- seurat@assays$RNA@data[genes, ]
  subset <- data.frame(data.matrix(t(subset)))
  seurat <- AddMetaData(seurat, metadata = subset)

  if(FALSE %in% gene.index){
      print("The following genes are not in the dataset: ")
      print(as.character(marker.genes$GENE)[!gene.index])
  }

########## Plot marker genes

  print("Plotting Marker Gene Expression...")
  
  meta <- seurat@meta.data
  sigs <- colnames(meta)[colnames(meta) %in% marker.genes$GENE]
  print(sigs)


    #### Tiff files of tSNE with marker gene expression in hexbin

    for (i in seq_along(sigs)) {

        goi <- sigs[i] 

        meta_ordered <- meta[order(meta[,goi], decreasing = F), ]

        plot.title <-  paste(goi, "\n", outName, "  ", nrow(meta), "cells")

        p <-  ggplot(meta_ordered, aes_string(x="tSNE_1", y="tSNE_2", z=goi)) +
                            stat_summary_hex(bins=100, fun = "mean") +
                             theme_classic(base_size=8) +
                             scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                             theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                            ggtitle(plot.title) + labs(z = goi)

        mypath <- paste0(markerGeneDir, outName, "_", goi, "_GeneMarkers_tSNE.tiff")

        tiff(file=mypath)
        print(p)
        dev.off()

    }


    #### Tiff files of umap with marker gene expression in hexbin

    for (i in seq_along(sigs)) {

        goi <- sigs[i] 

        meta_ordered <- meta[order(meta[,goi], decreasing = F), ]

        plot.title <-  paste(goi, "\n", outName, "  ", nrow(meta), "cells")

        p <-  ggplot(meta_ordered, aes_string(x="UMAP_1", y="UMAP_2", z=goi)) +
                            stat_summary_hex(bins=100, fun = "mean") +
                             theme_classic(base_size=8) +
                             scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                             theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                            ggtitle(plot.title) + labs(z = goi)

        mypath <- paste0(markerGeneDir, outName, "_", goi, "_GeneMarkers_UMAP.tiff")

        tiff(file=mypath)
        print(p)
        dev.off()

      }
}


StopWatchEnd$SaveNewMarkers <- Sys.time()
  


#########################################
# Subset cohort object into lineage-specific objects 
#########################################

print("")
print("*********************")
print("Subset and save lineages")
print(Sys.time())
print("********************")
print("")
StopWatchStart$SubsetSave <- Sys.time()

for (i in unique.lineages) {
  
    subset <- subset(seurat, subset = Subset_lineage == i)
    
    print(paste0("Total cells for lineage ", i, ": ", ncol(subset)))
    print("")
    print("*********************")
    print(paste0("Distribution of ", i, " across samples"))
    print(table(subset@meta.data$Sample_ID))
    print("********************")
    print("")

    dir.create(i)
  
    print(paste0("Saving Raw Counts for Lineage ", i, "...."))
  
    raw.count.dir <- paste0("./", i, "/", outName, "_", i, "_rawCounts")
    raw <- Matrix(subset@assays$RNA@counts, sparse = TRUE)
    write10xCounts(path = raw.count.dir,
                   x = raw
                   )
  
    print(paste0("Saving Normalized Counts for Lineage ", i, "...."))
    norm.count.dir <- paste0("./", i, "/", outName, "_", i, "_normCounts")  
    norm.count <- Matrix(subset@assays$RNA@data, sparse = TRUE)
    write10xCounts(path = norm.count.dir,
                   x = norm.count
                   )

    print(paste0("Saving Meta Data for Lineage ", i, "...."))
    meta <- data.frame(subset@meta.data)
    meta.file <- paste0("./", i, "/", outName, "_", i, "_metaData.csv")  
    write.csv(meta, file = meta.file)
    meta.file.2 <- paste0("./", i, "/", outName, "_", i, "_metaData.rds")  
    saveRDS(meta, file = meta.file.2)


    print(paste0("Saving Seurat Object for Lineage ", i, "...."))
    seurat.file <- paste0("./", i, "/", outName, "_", i, "_seurat.rds")  
    saveRDS(subset, file = seurat.file)

}


StopWatchEnd$SubsetSave <- Sys.time()




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
