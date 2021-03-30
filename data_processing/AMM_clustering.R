
###########
### Script for sample-level multi-resolution clustering of malignant cells from VkMYC mice with active-MM
### D.Croucher
### September 8, 2020
############


minResolution <- 0.4
maxResolution <- 1.4
fileName <- "MM1" #iterate for each sample
fdrCutoff <- 0.05
deGeneCutoff <- 5

#-----full seurat object 
subset <- readRDS("./path/to/mPC/seurat/object.rds")
subset <- subset(subset, subset = orig.ident == fileName)
dim(subset)


#-----run seurat pipeline 

#remove gene.filtering % of lowest library size

nCells <- mean((table(subset@meta.data$Sample_ID))) * 0.001
print(paste0("Cutoff: ", 0.001, " of average cell size (", round(nCells, 2), " cells)"))

seurat.obj <- CreateSeuratObject(counts = subset@assays$RNA@counts,
                                 meta.data = subset@meta.data[,c(1:19)],
                                 min.cells = nCells
                                )


seurat.obj <- NormalizeData(object = seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000
                           )

seurat.obj <- FindVariableFeatures(seurat.obj,
                                       selection.method = "vst",
                                       nfeatures = 3000
                                      )

seurat.obj <- ScaleData(seurat.obj,
                    vars.to.regress = "percent.mt",
                    features =  VariableFeatures(object = seurat.obj)
                   )


seurat.obj <- RunPCA(seurat.obj,
                     features = VariableFeatures(object = seurat.obj),
                     npcs = 100,
                     verbose = FALSE
                    )


pca <- data.frame(seurat.obj@reductions$pca@stdev)
pca$PC <- seq(1:nrow(pca))
pca <- pca[ ,c(2,1)]
colnames(pca)[2] <- "st.dev"
eigs <- pca$st.dev**2
pca$prop.var <- eigs / sum(eigs)

cutoff.point <- findCutoff(pca$PC[1:75],
                               pca$st.dev[1:75],
                               method="first",
                               0.01 #derivative cutoff
                              )
numPC <- round(cutoff.point$x)
pc.use <- numPC
pc.use

seurat.obj <- RunTSNE(seurat.obj, dims = 1:pc.use, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:pc.use, verbose = FALSE)

seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$umap@cell.embeddings)
                          )
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$tsne@cell.embeddings)
                          )

seurat.obj <- FindNeighbors(seurat.obj,
                            dims = 1:pc.use,
                            verbose = TRUE
                           )
DimPlot(seurat.obj)

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
                          )

##fix naming in meta.data

colnames(seurat.obj@meta.data) <- gsub("RNA_snn", "Seurat_cluster", colnames(seurat.obj@meta.data))
seurat.obj@meta.data$seurat_clusters <- NULL


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
markers <- do.call("rbind", markers)
rownames(markers) <- NULL

#remove genes that are above FDR cutoff
print(paste0("Filtering DE genes....FDR=", fdrCutoff))
markers_fdr <- markers[markers$p_val_adj <= fdrCutoff, ]
print(paste0("Removed ", nrow(markers) - nrow(markers_fdr), " genes...."))


# Sil Width
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

table(Idents(seurat.obj))
DimPlot(seurat.obj)
dim(seurat.obj)

print("Saving Seurat Object....")
seurat.file <- paste0("./data/", fileName, "_seurat.rds")
saveRDS(seurat.obj, file = seurat.file)




