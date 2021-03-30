# Figure 4

# Determine cluster-specific pathway terms using reactome (https://reactome.org/)
# Load reactome data and plot
```{r}


###########################################################
### Load DE results for each sample
###########################################################

reactome.files <- grep("reactome", list.files("./data/reactome"), value = T)

for (i in seq_along(reactome.files)) {
  
  if (i == as.integer(1)) {
    
    file <- paste0("./data/reactome/", reactome.files[i])
    gs <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    clustercolumn <- gsub("reactome_", "", reactome.files[i])
    clustercolumn <- gsub(".csv", "", clustercolumn)
    gs$ClusterID <- clustercolumn
    
    }
  
    if (i > as.integer(1)){
      
    file <- paste0("./data/reactome/", reactome.files[i])
    tmp <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    clustercolumn <- gsub("reactome_", "", reactome.files[i])
    clustercolumn <- gsub(".csv", "", clustercolumn)
    tmp$ClusterID <- clustercolumn
    
    gs <- rbind(gs, tmp)
    
    }
}

dim(gs)

length(unique(gs$ClusterID))

#remove entries with FDR > 0.05
gs <- gs[gs$Entities.FDR <= 0.05, ] #removes 4 cluster ids 

dim(gs)
length(unique(gs$ClusterID))

#remove entries with only 1 gene (Submitted.entities.found)
gs <- gs[grep(";", gs$Submitted.entities.found), ] 

library(stringr)
tmp <- data.frame(gs$Submitted.entities.found)
gs$count <- apply(tmp, 1, function(x) str_count(x, pattern = ";")) + 1
rm(tmp)

dim(gs)
length(unique(gs$ClusterID)) #removes an additional 1 cluster id 

###########################################################
### Create binary matrix of marker genes across clusters
###########################################################

b <- table(gs$Pathway.name, gs$ClusterID)

b <- as.matrix(b)
b[1:5, 1:5]

###########################################################
### Compute similarity and plot
###########################################################

#compute the overlap (similarity) of marker gene profiles for each cluster
#use the bianry distance in R, which is the same as computing the jaccard distance
#Jaccaard distance measures dissimilairty between sample sets, subtracting Jaccard coefficient from 1
# Jaccard distance = 1-(jaccard coefficient)

jaccard <- dist(t(b), method = "binary") 
str(jaccard)

hc <- hclust(jaccard, method = "complete")
plot(hc, hang = -1, cex = 0.4)

jaccard_index <- 1-as.matrix(jaccard)

library(ComplexHeatmap)
library(circlize)
library(viridis)



###########################################################
### Complex Heatmap
###########################################################


jaccard_index[1:5,1:5]


# set variables you need to add in annotations

goi <- c("Mki67", "Top2a")
sample <- c("MM1", "MM2", "MM3", "MM4", "MM5", "MM6", "MM7")

# 
runSignature <- "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY_MOUSE"
geneSignatures <- "~/OneDrive - UHN/Documents/10x_Experiments/MASTER_genesets.csv"
sigs <- read.csv(geneSignatures, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
runSig <- sigs[2:nrow(sigs), ,drop = F]
runSig <- runSig[, runSignature]
runSig <- runSig[!is.na(runSig)]
sigs <- list(runSig)
names(sigs) <- runSignature


meta <- read.csv("./data/scVKMYC_cohort_multires_harmony_theta05_AUCell_meta.csv", row.names = 1)
meta <- data.frame(meta$cellColor, 
                   row.names = rownames(meta))
nBreaks <- 5 # 
colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)



anno_data_full <- list()

for (i in sample) {
  
  file.name <- paste0("./data/", i, "_seurat.rds")
  
  dat <- readRDS(file.name)
  optimal.res <- paste0("Seurat_cluster_res.", unique(dat@meta.data$Optimal_res))
  dat@meta.data$cluster <- dat@meta.data[, optimal.res]
  
  
  # Add cell cycle genes
  
  dat <- ScaleData(dat,
                   vars.to.regress = "percent.mt",
                   features = goi)
  avg.scaled.exp <- data.frame(AverageExpression(dat, features = goi, slot = "scale.data"))
  colnames(avg.scaled.exp) <- paste0(i, "_c", levels(dat@meta.data$cluster))
  
  # Add cell cycle phase
  
  phase <- as.data.frame.matrix(prop.table(table(dat@meta.data$Phase, dat@meta.data$cluster), margin = 2))
  colnames(phase) <- paste0(i, "_c", levels(dat@meta.data$cluster))
  
  
  # Add chr5 status
  
  dat <- AddMetaData(dat, meta[rownames(meta) %in% rownames(dat@meta.data),], 
                     col.name = "cellColor")
  
  dat@meta.data[dat@meta.data$cellColor %in% colorPal_Pos[1:5], "chr5status"] <- "wt"
  dat@meta.data[dat@meta.data$cellColor %in% colorPal_Neg[1:5], "chr5status"] <- "del"
  
  chr5 <- as.data.frame.matrix(prop.table(table(dat@meta.data$chr5status, dat@meta.data$cluster), margin = 2))
  colnames(chr5) <- paste0(i, "_c", levels(dat@meta.data$cluster))

  
  
  # Add Protein Metabolism Signature score
  
  dat <- AddModuleScore(dat,
                        features = sigs,
                        ctrl = 25,
                        name = make.names(names(sigs))
  )
  start <- ncol(dat@meta.data)-length(sigs) + 1
  end <- ncol(dat@meta.data)
  colnames(dat@meta.data)[start:end] <- paste0(names(sigs), "_ModuleScore")
  dat@meta.data$module_sig <- dat@meta.data[,paste0(names(sigs), "_ModuleScore")]
  
  prot_meta_score <- dat@meta.data %>% dplyr::group_by(cluster) %>% 
    dplyr::select(module_sig) %>% dplyr::summarize_all(mean)
  prot_meta_score <- data.frame(prot_meta_score,
                                row.names = paste0(i, "_c", levels(dat@meta.data$cluster)))
  prot_meta_score$cluster <- NULL
  prot_meta_score <- t(prot_meta_score)
  
  anno_data <- rbind(avg.scaled.exp, phase, prot_meta_score, chr5)
  
  if (i == "MM1") {
    
    anno_data_full <- anno_data
    
  }else{
    
    anno_data_full <- cbind(anno_data_full, anno_data)
    
  }
  
}


#for paper - figure out average % in G2/M for clusters in SP2

reactome.protein.1 <- c("MM6_C1", "MM1_C8", "MM7_C2", "MM5_C1", "MM2_C2", "MM3_C5", "MM4_C0", "MM4_C7",
                        "MM4_C1", "MM1_C11", "MM7_C0") #0.6412949

reactome.prolif.2 <- c("MM4_C5", "MM4_C6", "MM5_C2", "MM6_C2", "MM4_C2", "MM2_C3", "MM7_C4", "MM1_C5", "MM3_C3") #0.6600436

exclude <- c(gsub("_C", "_c", reactome.protein.1), gsub("_C", "_c", reactome.prolif.2))

others <- colnames(anno_data_full)[!colnames(anno_data_full) %in% exclude] #0.05422651
reactome.prolif.2 <- gsub("_C", "_c", reactome.prolif.2)
reactome.protein.1 <- gsub("_C", "_c", reactome.protein.1)

reactome.protein.1
reactome.prolif.2
exclude

tmp1 <- t(anno_data_full)
mean(tmp1[reactome.protein.1, "G2M"]) #0.047314
sd(tmp1[reactome.protein.1, "G2M"]) #0.05649808

mean(tmp1[reactome.prolif.2, "G2M"]) #0.7399195
sd(tmp1[reactome.prolif.2, "G2M"]) #0.1116525

mean(tmp1[others, "G2M"]) #0.02194094
sd(tmp1[others, "G2M"]) #0.03412214

mean(tmp1[c(others, reactome.protein.1), "G2M"]) #0.02948428
sd(tmp1[c(others, reactome.protein.1), "G2M"]) #0.04281896

#--- resume plotting

#subset anno_data_full to include clusters in jaccard_index and reorder so it matches jaccard index

anno_data_full <- anno_data_full[, colnames(anno_data_full) %in% colnames(jaccard_index)]

jac_order <- colnames(jaccard_index)
tmp <- data.frame(t(anno_data_full))
tmp$Row.names <- rownames(tmp)
tmp$Row.names <- factor(tmp$Row.names,
                        levels = rownames(jaccard_index))
tmp <- tmp %>% arrange(Row.names)
rownames(tmp) <- tmp$Row.names

tmp$Row.names <- NULL
anno_data_full <- as.data.frame.matrix(tmp)
anno_data_full$Sample_ID <- gsub("_.*", "", rownames(anno_data_full))
anno_data_full$module_sig_scaled <- scale(anno_data_full$module_sig)


#select colours 
gg_color_hue <- function(n) { #Nina's function
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col <- c("#CFCFCF", "#909090", "#454545")

# For continuous values, the color mapping should be a color mapping function generated by circlize::colorRamp2()
library(circlize)
col_fun = colorRamp2(c(-2, 0, 5), c("blue", "white", "red"))


sample_colour_options <- brewer.pal(7, "Paired")

sample_colours <- list(Sample = c("MM1" = "#3A5FF6",
                                  "MM2" = "#FF56D1", #sub in phase/ch5 colours for legend 
                                  "MM3" = sample_colour_options[3],
                                  "MM4" = sample_colour_options[4],
                                  "MM5" = sample_colour_options[5],
                                  "MM6" = sample_colour_options[6],
                                  "MM7" = sample_colour_options[7]
)) #sub in phase/ch5 colours for legend 

chr5_col <- c("#3A5FF6", "#FF56D1")



#create annotations object
column_anno <- HeatmapAnnotation(Sample = anno_data_full[, "Sample_ID"], 
                                 Phase = anno_barplot(anno_data_full[, c("G1", "G2M", "S")],
                                                      gp = gpar(fill = col, col = "black")),
                                 Mki67 = anno_data_full[, "Mki67"],
                                 Top2a = anno_data_full[, "Top2a"],
                                 chr5 = anno_barplot(anno_data_full[, c("del", "wt")],
                                                      gp = gpar(fill = chr5_col, col = "black")),
                                 #ModuleScore = anno_data_full[, "module_sig"],
                                 gp = gpar(col = "black"),
                                 which = "column",
                                 gap = unit(2, "mm"),
                                 col = list(Sample = c("MM1" = sample_colour_options[1], #sub in phase/ch5 colours for legend 
                                                       "MM2" = sample_colour_options[2],
                                                       "MM3" = sample_colour_options[3],
                                                       "MM4" = sample_colour_options[4],
                                                       "MM5" = sample_colour_options[5],
                                                       "MM6" = sample_colour_options[6],
                                                       "MM7" = sample_colour_options[7]),
                                            Mki67 = col_fun,
                                            Top2a = col_fun,
                                            chr5 = c("DEL" = "#3A5FF6", 
                                                                 "WT" = "#FF56D1"))
)

column_reorder <- c(reactome.prolif.2, reactome.protein.1, others)
column_reorder <- column_reorder[column_reorder %in% colnames(jaccard_index)]
colnames(jaccard_index)

Heatmap(as.matrix(jaccard_index), 
        name = "Jaccard \nIndex",  #change this to "Phase" so you can get legend title
        col = c("white", magma(n=15, direction = -1)),
        top_annotation = column_anno,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        column_split = 4,
        row_split = 4,
        #column_order = column_reorder,
        #row_order = column_reorder
        
)




# add similiary-cluster relationship info to anno_data_full

anno_data_full$sample_cluster <- rownames(anno_data_full)

reactome.protein.1 <- c("MM6_C1", "MM1_C8", "MM7_C2", "MM5_C1", "MM2_C2", "MM3_C5", "MM4_C0", "MM4_C7",
                        "MM4_C1", "MM1_C11", "MM7_C0")
reactome.protein.1 <- gsub("_C", "_c", reactome.protein.1)

reactome.prolif.2 <- c("MM4_C5", "MM4_C6", "MM5_C2", "MM6_C2", "MM4_C2", "MM2_C3", "MM7_C4", "MM1_C5", "MM3_C3") #0.6600436
reactome.prolif.2 <- gsub("_C", "_c", reactome.prolif.2)

exclude <- c(reactome.protein.1, reactome.prolif.2)
others <- unique(gs$ClusterID)[!unique(gs$ClusterID) %in% exclude] 

anno_data_full[anno_data_full$sample_cluster %in% reactome.protein.1, "similarity_group"] <- "ProteinMetabolism"
anno_data_full[anno_data_full$sample_cluster %in% reactome.prolif.2, "similarity_group"] <- "CellCycle"
anno_data_full[anno_data_full$sample_cluster %in% others, "similarity_group"] <- "Dissimilar"

anno_data_full$similarity_group <- factor(as.character(anno_data_full$similarity_group), 
                                          levels = c("ProteinMetabolism", "CellCycle", "Dissimilar"))

d <- anno_data_full %>% group_by(similarity_group) %>% 
  select(module_sig) %>% summarize_all(mean)

ggplot(d, aes(x=similarity_group, y=module_sig)) + geom_point(fill = "grey", col = "black") +
  labs(y = "Mean ProteinMetabolism Signature Score") + theme_classic()
  
labels <- c("GCN2-ISR", "Cell Cycle", "Dissimilar")

ggplot(anno_data_full, aes(x=similarity_group, y=module_sig)) + geom_boxplot(fill = "grey", col = "black", outlier.shape = NA) +
  labs(y = "GCN2-ISR Signature Score", x= "Similarity Group") + theme_classic() +
  scale_x_discrete(labels= labels)
  











```
