###########################################################
#                Annotate cells with SingleR              #
#             D.Croucher (L.Richards template)            #
#                        Sept 2020                        #
###########################################################
### References:
### http://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
###

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load in seurat object + convert into sce object
### 2) Load in reference sets
### 3) Run SingleR cell annotation
### 4) Save results + plot
###########################################################

###########################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 48:00:00
## #SBATCH --mem=300G
## #SBATCH -p superhimem
## #SBATCH -c 50
## #SBATCH -N 1
## #SBATCH --account=pughlab
##
## module load R/3.6.1
##
## Rscript /cluster/projects/pughlab/projects/scVKMYC/scripts/singler_cellannotation.r 
## --seuratObj /cluster/projects/pughlab/projects/scVKMYC/analysis/test/seurat/scVKMYC_test_multires/data/scVKMYC_test_multires_seurat.rds \
## --outName scVKMYC_test \
## --referenceSets /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/singler/references/ImmGenData.rds \
## --referenceType bulkRNA \
## --logNormalize TRUE \
## --labelName label.main \
## --outputDir /cluster/projects/pughlab/projects/scVKMYC/analysis/test/singler/singlecell_mode \
## --scoreMode single

###########################################################

suppressMessages(library(optparse))
suppressMessages(library(SingleR))
suppressMessages(library(Seurat))
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(hash))
suppressMessages(library(rlist))

### custom function to load and rename Rdata obj
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


### custom function for running SingleR with high-dimensional datasets
runSinglerHD <- function(seurat_obj_in) {
  
  scData=seurat_obj_in; #seurat object
  predList <- list()
  
  print("")    
  print("***********************")
  print("Running modified Singler function to handle large datasets")
  cellsToPredict <- paste0(ncol(scData), " cells to predict")
  print(cellsToPredict)
  print(Sys.time())  
  print("***********************")  
    
    
  #Set up your intervals of 10,000 depending on how many cells you have
  
  splitNum <- as.integer(floor(ncol(scData)/10000)) #e.g. 20058 cells, splitNum = 2
  interval <- seq(0, splitNum*10000, by=10000) #[1]     0 10000 20000
  subsetIntervals <- list() 
    
  for (i in seq_along(1:(splitNum+1))) { #1 2 3
      
    if (!i == splitNum+1) { #first and second round
         
      subsetIntervals[[i]] <- noquote(paste(interval[i]+1, interval[i+1], sep = ":")) #1:10000 ... then 10001:20000      
        
    } else { #third round 
        
      subsetIntervals[[i]] <- noquote(paste((splitNum*10000+1), ncol(scData), sep = ":")) #20001:20058
        
    }


  #Subset seurat object and convert to sce

    cellSubset <- colnames(x = scData)[strsplit(as.character(subsetIntervals), 
                                     ":")[[i]][1]:strsplit(as.character(subsetIntervals), 
                                                           ":")[[i]][2]]
    seuratSubset <- subset(scData, cells =  cellSubset)
   
    sceSubset <- as.SingleCellExperiment(seuratSubset)
    sceSubset <- logNormCounts(sceSubset) #this should overwrite the data slot

  
  #Run singler on scDataSubset
  
  print("")    
  print("***********************")
  processName <- paste("Starting SingleR on cells", subsetIntervals[[i]])
  print(processName)
  print(Sys.time())  
  print("***********************")
    
  pred <- SingleR(test = sceSubset,
                  ref = ref,
                  labels = ref@colData[ ,labelName],
                  de.method = "classic",
                  check.missing = TRUE,
                  #clusters = clusters,
                  method = scoreMode)
   
     
  predList[[i]] <- pred
   
  print("")    
  print("***********************")
  processName <- paste("SingleR complete for cells", subsetIntervals[[i]])
  print(processName)
  print(Sys.time())  
  print("***********************")   
    
  }
  
  #Combine all pred dataframes and output
  
  print("")    
  print("***********************")
  processName <- paste("Combining SingleR predictions")
  print(processName)
  print(Sys.time())  
  print("***********************")   
    
  predFull = list.rbind(predList) 
    
  obj_out = predFull
                    
  return(obj_out) 
}



###########################################################
### USER DEFINED VARIABLES & PARSE OPTIONS
###########################################################

option_list <- list(make_option("--seuratObj",
                                type = "character",
                                default = NULL,
                                help = "path to seurat object to use as test dataset",
                                metavar= "character"
                               ),
                     make_option("--outName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files",
                                metavar= "character"
                               ),
                     make_option("--referenceSets",
                                type = "character",
                                help = "path to dir reference dataset, must be sce or se object and end in .rds. If multiple referneces are used, this should be a list object with list$References and list$Labels",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--referenceType",
                                type = "character",
                                help = "bulkRNA or scRNA; custom markerGenes not supported yet",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--outputDir",
                                type = "character",
                                help = "path to dir where all pipeline outputs will go",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--logNormalize",
                                type = "logical",
                                help = "TRUE/FALSE to run scater  LogNormCounts() on data",
                                default = FALSE,
                                metavar= "logical"
                                ),
                   make_option("--labelName",
                              type = "character",
                              help = "label metadata names in reference set for single reference sets (ie. sce$labelName)",
                              default = NULL,
                              metavar= "character"
                              ),
                   make_option("--scoreMode",
                                type = "character",
                                help = "can be 'single' or 'cluster';
                                right now only single mode supported",
                                default = FALSE,
                                metavar= "character"
                              ),
                    make_option("--clusterHeader",
                                type = "character",
                                help = "name of column in meta/colData to use if scoring in cluster mode",
                                default = FALSE,
                                metavar= "character"
                              )
                  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

seuratObj <- opt$seuratObj
outName <- opt$outName
referenceSets <- opt$referenceSets
referenceType <- opt$referenceType
logNormalize <- opt$logNormalize
labelName <- opt$labelName
scoreMode <- opt$scoreMode
clusterHeader <- opt$clusterHeader
outputDir <- opt$outputDir



###########################################################
### 1) Set up input data
###########################################################

### set working dir
setwd(outputDir)
dir.create(outName)
setwd(outName)


print("")
print("***********************")
print("Loading data")
print(Sys.time())
print("***********************")

print("Loading seurat object...")
##read in data
dat <- readRDS(seuratObj)

if (class(dat) == "Seurat"){

    print("Converting to SingleCellExperiment")
    #print(paste0("Extracting slot... ", seuratExtract))
    sce <- as.SingleCellExperiment(dat)

} else { print("Error: Input file not Seurat Object") }

print(paste0("Dimensions of input data...",
              nrow(dat),
              " genes x ",
              ncol(dat),
              " cells"
            )
    )



###normalize data if needed
if(logNormalize == TRUE){

  print("Running LogNormCounts()...")
  #sce <- SingleCellExperiment(assays=list(counts = dat))
  #this is automatically performed on count matrix
  sce <- logNormCounts(sce) #this should overwrite the data slot

} else if (logNormalize == FALSE){

  print("Using raw or pre-normalized counts...")

}

##########################################################
### 2) Set up reference data
##########################################################

print("")
print("***********************")
print("Setting up reference data")
print(Sys.time())
print("***********************")

print("Loading reference set(s)....")
ref <- readRDS(referenceSets)

if(class(ref) == "list"){

    print(paste0("Multiple references loaded....n=", length(ref$References)))
    print(names(ref$References))
    print(cat(unique(unlist(ref$Labels)), sep = "\n"))

} else {

  print(paste0("One reference loaded... ", basename(referenceSets)))
  print(paste0("Labels being used....n=", length(unique(ref@colData[ ,labelName]))))
  print(cat(unique(ref@colData[ ,labelName]), sep = "\n"))

}


###########################################################
### 3) Run SingleR annotation
###########################################################

print("")
print("***********************")
print("Running SingleR HD")
print(Sys.time())
print("***********************")

#coldata <- data.frame(colData(sce))
#coldata$cluster <- coldata[ ,colnames(coldata) == paste0("Seurat_cluster_res.", unique(coldata$Optimal_res))]
#clusters <- factor(coldata[ ,clusterHeader])

#only set up to run with bulk data using single reference

print("Bulk RNA-seq reference data set(s)....")

print(paste0("Running in ", scoreMode, " mode with one reference...."))
    
pred <- runSinglerHD(dat)



###########################################################
### 4) Save SingleR results
###########################################################

print("")
print("***********************")
print("Saving SingleR prediction results")
print(Sys.time())
print("***********************")

### 5.1) Save SingleR output
ref.name <- gsub(".rds", "", basename(referenceSets))
save.name <- paste0(outName, ".", ref.name, ".SingleR_pred.rds")
print(save.name)
saveRDS(pred, file = save.name)



###########################################################
### 5) Plot data onto UMAP / tSNE
###########################################################

print("")
print("***********************")
print("Plot Results")
print(Sys.time())
print("***********************")

###########################
### 5.1) Merge SingleR results with Seurat metadata
###########################
print("Merging SingleR labels and scores with Seurat metadata...")
### add labels

if (scoreMode == "single"){

  dat@meta.data$SingleR_FirstLabels <- pred$first.labels
  dat@meta.data$SingleR_PrunedLabels <- pred$pruned.labels
  dat@meta.data$SingleR_Labels <- pred$labels
  print("Agreement between final and pruned labels....")
  table(dat@meta.data$SingleR_Labels == dat@meta.data$SingleR_PrunedLabels)
  ### add SingleR scores
  scores <- data.frame(pred$scores)
  rownames(scores) <- rownames(dat@meta.data)
  colnames(scores) <- paste0("SingleR_Score.", colnames(scores))
  dat <- AddMetaData(dat,
                   metadata = scores
                 )

} else if (scoreMode == "cluster"){
  #extract optimal clustering
  opt_clust_col <- paste0("Seurat_cluster_res.",
                          unique(dat@meta.data$Optimal_res))
  #add C and +1 to cluster number to be compatible with plotting
  #factor levels into right order and reset ident
  dat[["SingleR_ClusterID"]] <- paste0("C",
                                      as.numeric(dat@meta.data[ ,opt_clust_col]))
  dat@meta.data$SingleR_ClusterID <- factor(dat@meta.data$SingleR_ClusterID,
                                        levels = paste0("C", 1:length(unique(dat@meta.data$SingleR_ClusterID))))

  # Single R first labels
  Idents(object = dat) <- "SingleR_ClusterID"
  first.labels <- pred$first.labels
  names(first.labels) <- levels(dat)
  dat <- RenameIdents(dat, first.labels)
  dat[["SingleR_FirstLabels"]] <- Idents(dat)

  # Single R first labels
  Idents(object = dat) <- "SingleR_ClusterID"
  pruned.labels <- pred$pruned.labels
  names(pruned.labels) <- levels(dat)
  dat <- RenameIdents(dat, pruned.labels)
  dat[["SingleR_PrunedLabels"]] <- Idents(dat)

  # Single R labels
  Idents(object = dat) <- "SingleR_ClusterID"
  labels <- pred$labels
  names(labels) <- levels(dat)
  dat <- RenameIdents(dat, labels)
  dat[["SingleR_Labels"]] <- Idents(dat)


}


###########################
### 5.2) Default SingleR plots
###########################

print("Plotting SingleR default plots.....")

if (ncol(dat) < 30000) {

    singleR.plots <- paste0(outName, ".SingleR_heatmap.pdf")
    pdf(singleR.plots, width = 20, height = 15)
    if (scoreMode == "single"){
      plotScoreHeatmap(pred,
                    annotation_col=as.data.frame(colData(sce)[ ,c("Sample_ID", "Cohort_ID", "Sex"),drop=FALSE]
                    ))
    } else if (scoreMode == "cluster"){
      plotScoreHeatmap(pred)
    }
    dev.off()

}else{
    

    #subsample so less than 30000 cells are included (evenly distributed across samples)
    print("Plotting SingleR default heatmap with 30,000 cells only.....")

        cells.per.sample <- 30000/length(unique(dat@meta.data$Sample_ID))
    
        Idents(dat) <- "Sample_ID" 
    
        counts <- data.frame(table(dat@meta.data$Sample_ID)) #how many cells in each sample
        rownames(counts) <- counts$Var1 #set rownames to sample_ID
        split.variable <- as.character(counts$Var1) #unique sample_IDs

        cells <- list()
                
            for (i in seq_along(split.variable)) {
                countsTmp <- counts[i, ]
                if (countsTmp$Freq > cells.per.sample) {
                    DS <- sample(1:length(WhichCells(dat, ident = split.variable[i])), cells.per.sample, replace = FALSE)
                    cells[[i]] <- WhichCells(dat, ident = split.variable[i])[DS]
                } else {
                    cells[[i]] <- WhichCells(dat, ident = split.variable[i])
                }
            }
    cells_unlist <- unlist(cells)
    
    #subset pred object
    predDS <- pred[rownames(pred) %in% cells_unlist, ]
    
    #subset sce object
    sceDS <- sce[, colnames(sce) %in% cells_unlist]
    
    #plot heatmap of downsampled data
    singleR.plots <- paste0(outName, ".SingleR_heatmap.pdf")
    pdf(singleR.plots, width = 20, height = 15)
    if (scoreMode == "single"){
      plotScoreHeatmap(predDS,
                    annotation_col=as.data.frame(colData(sceDS)[ ,c("Sample_ID", "Cohort_ID", "ident"),drop=FALSE]
                    ))
    } else if (scoreMode == "cluster"){
      plotScoreHeatmap(predDS)
    }
    dev.off()
    
}

  
    
singleR.plots2 <- paste0(outName, ".SingleR_pruneViolins.pdf")
pdf(singleR.plots2, width = 20, height = 15)
plotScoreDistribution(pred,
                      show = "delta.med",
                      ncol = 3,
                      show.nmads = 3
                    )
dev.off()


###########################
### 5.3) Average cluster Heatmap
###########################

meta <- dat@meta.data
meta$cluster <- meta[ ,colnames(meta) == paste0("Seurat_cluster_res.", unique(meta$Optimal_res))]
meta$cluster <- factor(as.character(meta$cluster),
                       levels = c(1:max(as.numeric(meta$cluster))))

if (scoreMode == "single"){
  print("Plotting cluster heatmap...")
  cent <- meta %>% group_by(cluster) %>% select(colnames(meta)[grep("SingleR_Score", colnames(meta))]) %>% summarize_all(mean)
  cent <- data.frame(cent)
  rownames(cent) <- cent$cluster
  cent$cluster <- NULL
  colnames(cent) <- gsub("SingleR_Score.", "", colnames(cent))
  cent <- scale(cent)
  cent[cent > 2] <- 2
  cent[cent < -2] <- -2
  heatmap.name <- paste0(outName, ".SingleR_ClusterHeatmap.pdf")
  pheatmap(t(cent), file = heatmap.name)
}

###########################
### 5.4) Plot proportion cell types
###########################

### 5.4.1) all labels
print("Plotting proportions...")
props <- prop.table(table(meta$SingleR_Labels, meta$cluster), margin = 2)
props <- data.frame(data.matrix(props))
colnames(props) <- c("CellType", "Cluster", "Proportion")
cols <- c(brewer.pal(n = 8, name = "Dark2"),
          brewer.pal(n = 12, name = "Paired"),
          brewer.pal(n = 8, name = "Set2")
         )
p4 <- ggplot() + geom_bar(aes(y = Proportion, x = Cluster, fill = CellType), data = props, stat="identity") +
      ggtitle("All Labels") + scale_fill_manual(values=cols)

### 5.4.2) collapse labels down

#collapse <- strsplit(as.character(meta$SingleR_Labels), "_", fixed = TRUE)
#meta$SingleR_CollapsedLabels <- sapply(collapse, "[[", 2)
#meta$SingleR_CollapsedLabels <- gsub("^Macrophage$", "Microglia/Macrophage", meta$SingleR_CollapsedLabels)
#meta$SingleR_CollapsedLabels <- gsub("^Microglia$", "Microglia/Macrophage", meta$SingleR_CollapsedLabels)
#meta$SingleR_CollapsedLabels <- gsub("^Myeloid$", "Microglia/Macrophage", meta$SingleR_CollapsedLabels)
#props <- prop.table(table(meta$SingleR_CollapsedLabels, meta$cluster), margin = 2)
#props <- data.frame(data.matrix(props))
#colnames(props) <- c("CellType", "Cluster", "Proportion")
#cols <- c(brewer.pal(n = 8, name = "Dark2"),
#          brewer.pal(n = 12, name = "Paired"),
#          brewer.pal(n = 8, name = "Set2")
#         )
#p5 <- ggplot() + geom_bar(aes(y = Proportion, x = Cluster, fill = CellType), data = props, stat="identity") +
#      ggtitle("Condensed Labels") + scale_fill_manual(values=cols)

prop.file <- paste0(outName, ".SingleR_CellTypeProportions.pdf")
pdf(prop.file, width = 10, height = 6)
p4
#p5
dev.off()

###########################
### 5.5) Assign each cluster to max proportion
###########################
if (scoreMode == "single"){
  max.label <- props %>% group_by(Cluster) %>%
  select(Proportion) %>% summarize_all(max)
  rownames(max.label)
  assign.labs <- props[props$Proportion %in% max.label$Proportion, ]
  assign.labs

  Label_MaxProp <- c()
  for (i in 1:nrow(meta)){
      clust <- meta[i, ]$cluster
      Label_MaxProp[i] <- as.character(assign.labs[assign.labs$Cluster == clust, ]$CellType)
}
  meta$Label_MaxProp <- Label_MaxProp
}

###########################
### 5.5) UMAPs colored by each signature
###########################
print("Plotting UMAPs....")
### 5.5.1) plot transcriptional cluster
cent <- meta %>% group_by(cluster) %>% select(UMAP_1,
    UMAP_2) %>% summarize_all(median)
umap.c <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "cluster")) +
          geom_point(alpha = 0.3) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs("Cluster") + ggtitle(paste0(outName, " - Cluster (", nrow(meta), " cells)")) +
          geom_label_repel(aes(label = cluster), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none')

### 5.5.2) plot by Sample_ID
cent <- meta %>% group_by(Sample_ID) %>% select(UMAP_1,
    UMAP_2) %>% summarize_all(median)
umap.p <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Sample_ID")) +
          geom_point(alpha = 0.3) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs("Cluster") + ggtitle(paste0(outName, " - Sample (", nrow(meta), " cells)")) +
          geom_label_repel(aes(label = Sample_ID), data = cent, show.legend = F, size = 2)  +
           theme(legend.position='none')

### 5.5.3) plot SingleR Collapsed labels

#cent <- meta %>% group_by(SingleR_CollapsedLabels) %>% select(UMAP_1,
#    UMAP_2) %>% summarize_all(median)
#umap.l <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "SingleR_CollapsedLabels")) +
#          geom_point(alpha = 0.3) + theme_classic() +
#          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#          labs("Cluster") + ggtitle(paste0(outName, " - SingleR Collapsed Labels (", nrow(meta), " cells)")) +
#          geom_label_repel(aes(label = SingleR_CollapsedLabels), data = cent, show.legend = F, size = 2)  +
#           theme(legend.position='none')

### 5.5.4) plot SingleR Assigned labels
if (scoreMode == "single"){

  cent <- meta %>% group_by(SingleR_Labels) %>% select(UMAP_1,
    UMAP_2) %>% summarize_all(median)
  umap.a <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "SingleR_Labels")) +
          geom_point(alpha = 0.3) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs("Cluster") + ggtitle(paste0(outName, " - SingleR Assigned Labels (", nrow(meta), " cells)")) +
          geom_label_repel(aes(label = SingleR_Labels), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none')

}

if(scoreMode == "single"){

  umap.file <- paste0(outName, ".SingleR_UMAPs.pdf")
  pdf(umap.file, width = 12, height = 12)
  print(ggarrange(umap.c, umap.p, umap.a, ncol = 2, nrow = 2))
  dev.off()

} else if (scoreMode == "cluster"){

  umap.file <- paste0(outName, ".SingleR_UMAPs.pdf")
  pdf(umap.file, width = 12, height = 12)
  print(ggarrange(umap.c, umap.p, ncol = 2, nrow = 2))
  dev.off()

}


###########################
### 5.6) tSNEs colored by each signature
###########################
print("Plotting tSNEs....")
### 5.5.1) plot transcriptional cluster
cent <- meta %>% group_by(cluster) %>% select(tSNE_1,
    tSNE_2) %>% summarize_all(median)
tsne.c <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "cluster")) +
          geom_point(alpha = 0.3) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs("Cluster") + ggtitle(paste0(outName, " - Cluster (", nrow(meta), " cells)")) +
          geom_label_repel(aes(label = cluster), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none')

### 5.5.2) plot by Sample_ID
cent <- meta %>% group_by(Sample_ID) %>% select(tSNE_1,
    tSNE_2) %>% summarize_all(median)
tsne.p <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Sample_ID")) +
          geom_point(alpha = 0.3) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs("Cluster") + ggtitle(paste0(outName, " - Sample (", nrow(meta), " cells)")) +
          geom_label_repel(aes(label = Sample_ID), data = cent, show.legend = F, size = 2)  +
           theme(legend.position='none')

### 5.5.3) plot SingleR Collapsed labels

#cent <- meta %>% group_by(SingleR_CollapsedLabels) %>% select(UMAP_1,
#    UMAP_2) %>% summarize_all(median)
#umap.l <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "SingleR_CollapsedLabels")) +
#          geom_point(alpha = 0.3) + theme_classic() +
#          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#          labs("Cluster") + ggtitle(paste0(outName, " - SingleR Collapsed Labels (", nrow(meta), " cells)")) +
#          geom_label_repel(aes(label = SingleR_CollapsedLabels), data = cent, show.legend = F, size = 2)  +
#           theme(legend.position='none')

### 5.5.4) plot SingleR Assigned labels
if (scoreMode == "single"){

  cent <- meta %>% group_by(SingleR_Labels) %>% select(tSNE_1,
    tSNE_2) %>% summarize_all(median)
  tsne.a <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "SingleR_Labels")) +
          geom_point(alpha = 0.3) + theme_classic() +
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
          labs("Cluster") + ggtitle(paste0(outName, " - SingleR Assigned Labels (", nrow(meta), " cells)")) +
          geom_label_repel(aes(label = SingleR_Labels), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none')

}

if(scoreMode == "single"){

  tsne.file <- paste0(outName, ".SingleR_tSNEs.pdf")
  pdf(tsne.file, width = 12, height = 12)
  print(ggarrange(tsne.c, tsne.p, tsne.a, ncol = 2, nrow = 2))
  dev.off()

} else if (scoreMode == "cluster"){

  tsne.file <- paste0(outName, ".SingleR_tSNEs.pdf")
  pdf(tsne.file, width = 12, height = 12)
  print(ggarrange(tsne.c, tsne.p, ncol = 2, nrow = 2))
  dev.off()

}
###########################################################
### 6) Save metadata results
###########################################################

print("")
print("***********************")
print("Saving meta data")
print(Sys.time())
print("***********************")
save.name2 <- paste0(outName, ".", ref.name, ".SingleR_metadata.rds")
print(save.name2)
saveRDS(meta, file = save.name2)

print("Saving Seurat Object....")
seurat.file <- paste0(outName, ".", ref.name, "_SinglerAnnotations_seurat.rds")
saveRDS(dat, file = seurat.file)

###########################################################
print("")
print("***********************")
print("END OF SCRIPT")
print(Sys.time())
print("***********************")
print("")
print(sessionInfo())
###########################################################
