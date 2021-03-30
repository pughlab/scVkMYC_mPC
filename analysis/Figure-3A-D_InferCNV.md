# Fig 3 iCNV
## Run InferCNV

```{r}

#locally from R

# https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
# https://github.com/broadinstitute/inferCNV/wiki/File-Definitions

library(infercnv)
library(Seurat)

setwd(outputDir)

#load full data set 
seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_PConly", "_seurat.rds")
dat <- readRDS(seurat.file)

#create a column for icnv_ident
dat@meta.data$malig_status <- Idents(dat)
dat@meta.data$icnv_ident <- paste0("malignant_", dat@meta.data$Sample_ID)
dat@meta.data[dat@meta.data$malig_status == "nPC", "icnv_ident"] <- "nPC"
Idents(dat) <- "icnv_ident"

#Raw Counts Matrix for Genes x Cells
counts.matrix <- as.matrix(dat@assays$RNA@counts)

counts.matrix.file.name <- paste0("./data/tsne_refined/icnv/", "scVKMYC_TSNErefined_PConly.matrix")
write.table(counts.matrix,
            file=counts.matrix.file.name,
            quote=F, 
            sep="\t")

#Sample annotation file
anno <- data.frame(Idents(dat), row.names = names(Idents(dat)))

anno.file.name <- paste0("./data/tsne_refined/icnv/", "scVKMYC_TSNErefined_PConly_anno.txt")

write.table(anno,
            file=anno.file.name,
            row.names = T,
            quote=F, sep= "\t", col.names = FALSE)


#Gene ordering file (created as per InferCNV wiki instructions)

gene.order.file.name <- "./data/tsne_refined/icnv/custom_2_mm10_hMYC_v2_gen_pos_nodup_SORTEDchr.txt"


#Run inferCNV

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts.matrix.file.name,
                                    annotations_file=anno.file.name,
                                    delim="\t",
                                    gene_order_file=gene.order.file.name,
                                    ref_group_names=c("nPC")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=0.1 works well for 10x Genomics
                             out_dir="./data/tsne_refined/icnv/scVKMYC_TSNErefined_PConly_HMMM", 
                             HMM_type = 'i3',
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE, 
                             num_threads = 4, #only 4 cores available 
                             analysis_mode='subclusters')

```

## Average CNV-subclones

```{r}

###########################################################
### create cell.groupings file for all samples 
###########################################################

icnv_output <- "./data/tsne_refined/icnv/scVKMYC_TSNErefined_PConly_HMMM/"

cell.groupings <- read.table(paste0(icnv_output, grep("cell_groupings", list.files(icnv_output), value = T)), header = T)

cell.groupings <- cell.groupings[!cell.groupings$cell_group_name %in% grep("nPC", cell.groupings$cell_group_name, value = T), ]

cell.groupings$sample <- gsub(".malignant.*", "", cell.groupings$cell_group_name)
cell.groupings$sample <- gsub("malignant_", "", cell.groupings$sample)

rownames(cell.groupings) <- cell.groupings$cell
head(cell.groupings)


cell.groupings.file <- "./data/tsne_refined/icnv/icnv_cell_groupings.txt"

write.table(cell.groupings, file = cell.groupings.file)


for (i in unique(cell.groupings$sample)) {
  
  if (i == unique(cell.groupings$sample)[1]) {
  
    cell.groupings.subset <- cell.groupings[cell.groupings$sample == i, ]
    
    for (k in seq_along(unique(cell.groupings.subset$cell_group_name))) {
    
      subclone.name <- as.character(unique(cell.groupings.subset$cell_group_name))[k]
      cell.groupings.subset[cell.groupings.subset$cell_group_name == subclone.name, "subclone"] <- k
      
    
    }
    
  final.cell.groupings <- cell.groupings.subset

  } else {
    
    cell.groupings.subset <- cell.groupings[cell.groupings$sample == i, ]
    
    for (k in seq_along(unique(cell.groupings.subset$cell_group_name))) {
    
    subclone.name <- as.character(unique(cell.groupings.subset$cell_group_name))[k]
    cell.groupings.subset[cell.groupings.subset$cell_group_name == subclone.name, "subclone"] <- k
    
    }
    
  final.cell.groupings <- rbind(final.cell.groupings, cell.groupings.subset)

  }
  
}

cell.groupings <- final.cell.groupings
cell.groupings$subclone <- as.character(cell.groupings$subclone)

cell.groupings.file <- "./data/tsne_refined/icnv/icnv_cell_groupings.txt"

write.table(cell.groupings, file = cell.groupings.file)

###########################################################
### Add to seurat metadata
###########################################################

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_mPConly", "_seurat.rds")
dat <- readRDS(seurat.file)

cell.groupings$cell <- NULL
cell.groupings$sample <- NULL


head(cell.groupings)

dat <- AddMetaData(dat, metadata = cell.groupings)

meta <- dat@meta.data


ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "subclone")) +
  geom_point(size = 1, alpha = 0.8) + theme_classic() +
  theme_bw() + theme(axis.text.x = element_blank(), 
                     axis.text.y = element_blank(), 
                     axis.ticks = element_blank(),
                     panel.border = element_rect(linetype = "solid",
                                                fill = NA, size = 1),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     legend.position = "none") +
  labs(x="UMAP_1", y="UMAP_2") +
  facet_grid(cluster~Sample_ID) +
  scale_color_brewer(palette = "Set1", 8)



```

# iCNV Complex heatmap
```{r}

#after running iCNV 

setwd(outputDir)
icnv.dir <- "./data/tsne_refined/icnv/scVKMYC_TSNErefined_PConly_HMMM/"

observations <- read.csv(paste0(icnv.dir, "infercnv.observations.txt"), sep = " ")
references <- read.csv(paste0(icnv.dir, "infercnv.references.txt"), sep = " ")

pos_file <- read.csv("./data/tsne_refined/icnv/custom_2_mm10_hMYC_v2_gen_pos_nodup_SORTEDchr.txt", 
                     sep = "\t",
                     header = FALSE)


###########################################################
### order by genome coordinates/sub-in chr number (using pos_file)
###########################################################


ordered_genes <- pos_file[pos_file$V1 %in% rownames(observations), ]

observations$Row.names <- rownames(observations)

observations$Row.names <- factor(observations$Row.names,
                                 levels = ordered_genes$V1)
observations <- observations %>% arrange(Row.names)
observations$Row.names <- NULL


observations[1:5,1:5]
head(ordered_genes)

ordered_genes$V2 <- factor(as.character(ordered_genes$V2),
                           levels = c(paste0("chr", rep(1:19) )))

#cells = rows, genes = columns
observations_mat <- t(as.matrix(observations))



###########################################################
### Bring in meta data
###########################################################

cell.groupings.file <- "./data/tsne_refined/icnv/icnv_cell_groupings.txt"

cell.groupings <- read.table(cell.groupings.file, header = T, row.names = 1)

head(cell.groupings)

#add subclone groupings to observations_mat (matching by cell id)
observations_mat[1:5,1:5]
dim(observations_mat)

observations_mat_subclones <- merge(observations_mat, cell.groupings, by = 0 , all = T)
rownames(observations_mat_subclones) <- observations_mat_subclones$Row.names
observations_mat_subclones$Row.names <- NULL
observations_mat_subclones[1:5,6218:6223]

subclone_mean_observations <- observations_mat_subclones
subclone_mean_observations$cell <- NULL
subclone_mean_observations$sample <- NULL
subclone_mean_observations$subclone <- NULL

#compute mean for each subclone
subclone_mean_observations <- subclone_mean_observations %>% group_by(cell_group_name) %>% 
  summarise_all(mean)


#reorder/regroup
subclone_mean_observations$sample <- gsub(".malignant.*", 
                                          "", 
                                          subclone_mean_observations$cell_group_name)
subclone_mean_observations$sample <- gsub("malignant_", "", subclone_mean_observations$sample)

subclone_mean_observations$sample <- factor(as.character(subclone_mean_observations$sample),
                         levels = c(
                                     "ND1",  "ND4", "ND5",
                                     "AMG1", "AMG2", "AMG3",
                                     "MM1", "MM2", "MM3","MM4","MM5","MM6","MM7")) 

subclone_mean_observations <- subclone_mean_observations %>% arrange(sample)
subclone_mean_observations <- data.frame(subclone_mean_observations)

rownames(subclone_mean_observations) <- subclone_mean_observations$cell_group_name
subclone_mean_observations$cell_group_name <- NULL


#use this to subset rows of heatmap
row_split <- subclone_mean_observations
row_split[1:5,6218:6220]

subclone_mean_observations$sample <- NULL



###########################################################
### Set up colours 
###########################################################

#Heatmap 
cols <- rev(c("#67001f",
          "#b2182b",
          "#d6604d",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          ))
colfunc <- colorRampPalette(cols)
heatmap.cols <- colfunc(50)

#Chromosomes 
chr.colours <- c(rep(c("lightgrey", "darkgrey"), 9), "lightgrey")




###########################################################
### Specify annotation colours 
###########################################################

anno_colours_column <- list(Chromosome = chr.colours)
names(anno_colours_column$Chromosome) <- levels(ordered_genes$V2)


###########################################################
### Create annotations object
###########################################################

library(ComplexHeatmap)

column_annotations <- HeatmapAnnotation(Chromosome = ordered_genes$V2,
                                        which = "column", 
                                        border = TRUE, 
                                        col = anno_colours_column,
                                        annotation_legend_param = list(
                                          Chromosome = list(nrow = 1)))


###########################################################
### Plot heatmap, apply clustering to rows (cells)
###########################################################


Heatmap(subclone_mean_observations, 
        name = "mat", 
        #top_annotation = column_annotations,
        #right_annotation = row_anotations,
        #show_row_dend = TRUE,
        show_column_dend = FALSE,
        border = TRUE, 
        row_names_gp = gpar(fontsize = 12),
        cluster_rows = FALSE,
        #clustering_distance_rows = "euclidean",
        #clustering_method_rows = "ward.D",
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        col = heatmap.cols,
        cluster_column_slices = FALSE,
        column_split = c(ordered_genes$V2),
        row_split = c(row_split$sample), 
        row_dend_width = unit(2, "cm")
        )







```


## Scale row height
```{r}

###########################################################
### Find subclone breakdown for each tumour
###########################################################

# % of cells per clone

df <- data.frame(prop.table(x = table(cell.groupings$cell_group_name, cell.groupings$sample),
                                       margin = 2))


df <- df[!df$Var2 == "ND3", ]
df[df==0] <- NA
df <- df[!is.na(df$Freq), ]

df$Var2 <- factor(as.character(df$Var2),
                                         levels = rev(c("ND1", "ND4", "ND5",
                                                    "AMG1", "AMG2", "AMG3",
                                                    "MM1", "MM2", "MM3","MM4","MM5","MM6","MM7")))
df <- df %>% arrange(desc(Freq)) %>% arrange(Var2) 


df[df$Var2 %in% grep("ND", df$Var2, value = T), "Cohort_ID"] <- "ND"
df[df$Var2 %in% grep("AMG", df$Var2, value = T), "Cohort_ID"] <- "AMG"
df[df$Var2 %in% grep("MM", df$Var2, value = T), "Cohort_ID"] <- "MM"
df$Cohort_ID <- factor(df$Cohort_ID, 
                       levels = c("ND", "AMG", "MM"))






###########################################################
### Integrate that with observation data
###########################################################

subset <- subclone_mean_observations

subset$clone <- rownames(subset)
subset[,1:5]


#Add that info to subset
df <- df[df$Var1 %in% rownames(subset), ]
rownames(df) <- df$Var1
df

subset <- merge(df, subset, by = 0, all = T)

#Clean up and re-name columns
subset$Row.names <- NULL
colnames(subset) <- c("clone", "sample", "clone_freq", c(colnames(subset)[4:ncol(subset)]))

melted_subset <- melt(subset, id.vars = c("clone", "sample", "clone_freq"))
head(melted_subset)



###########################################################
### Repeat each clone by their proportion
###########################################################

#round proportion first (use 1000 because you lose clones by rounding down if you use 100)

subset_rounded <- subset
subset_rounded$clone_freq <- subset_rounded$clone_freq*1000

subset_rounded$clone_freq <- round(subset_rounded$clone_freq, digits = 0)
subset_rounded[,1:5]

subset_rounded_repeated <- subset_rounded

for (i in 1:nrow(subset_rounded)) {
  
  subset_rounded_repeated_tmp <- do.call("rbind", 
                                     replicate(subset_rounded[i,3], 
                                               subset_rounded[i, ], 
                                               simplify = FALSE))

  subset_rounded_repeated_tmp <- data.frame(subset_rounded_repeated_tmp)
  
  subset_rounded_repeated_tmp[1:5,1:5]
 
  if (i == as.integer(1)){

      subset_rounded_repeated <- subset_rounded_repeated_tmp
    
    }
  
    if (i > as.integer(1)){
      
      subset_rounded_repeated <- rbind(subset_rounded_repeated, subset_rounded_repeated_tmp)
      subset_rounded_repeated[, 1:5]
    }
  
  }
subset_rounded_repeated_stash <- subset_rounded_repeated

#arrange by sample (ND at top)

subset_rounded_repeated$sample <- factor(as.character(subset_rounded_repeated$sample),
                         levels = c( "ND1",  "ND4", "ND5",
                                     "AMG1", "AMG2", "AMG3",
                                     "MM1", "MM2", "MM3","MM4","MM5","MM6","MM7")) 

subset_rounded_repeated <- subset_rounded_repeated %>% arrange(sample)
subset_rounded_repeated <- data.frame(subset_rounded_repeated)

#create row split for heatmap 

row_split <- subset_rounded_repeated
row_split[1:5,1:5]


#finalize heatmap matrix

subset_rounded_repeated$sample <- NULL
subset_rounded_repeated$clone_freq <- NULL
subset_rounded_repeated$clone <- NULL
subset_rounded_repeated$clone.1 <- NULL


#Heatmap colours 


#Chromosomes 
chr.colours <- c(rep(c("lightgrey", "darkgrey"), 9), "lightgrey")



###########################################################
### Specify annotation colours 
###########################################################

anno_colours_column <- list(Chromosome = chr.colours)
names(anno_colours_column$Chromosome) <- levels(ordered_genes$V2)


###########################################################
### Create annotations object
###########################################################

library(ComplexHeatmap)

column_annotations <- HeatmapAnnotation(Chromosome = ordered_genes$V2,
                                        which = "column", 
                                        border = TRUE, 
                                        col = anno_colours_column,
                                        annotation_legend_param = list(
                                          Chromosome = list(nrow = 1)))



Heatmap(as.matrix(subset_rounded_repeated), 
        name = "mat", 
        top_annotation = column_annotations,
        #left_annotation = row_annotations,
        #right_annotation = row_anotations,
        #show_row_dend = TRUE,
        #show_column_dend = FALSE,
        border = TRUE, 
        row_names_gp = gpar(fontsize = 8),
        cluster_rows = FALSE,
        #clustering_distance_rows = "euclidean",
        #clustering_method_rows = "ward.D",
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        col = heatmap.cols,
        cluster_column_slices = FALSE,
        column_split = c(ordered_genes$V2),
        row_split = c(row_split$sample), 
        #row_dend_width = unit(2, "cm")
        )





```


# Map CNVs to transcriptional clusters
```{r}

#need to run sample-specific clustering first


#Tumours are comprised of subpopulations with distinct CNV profiles
#Do tumours also demonstrate transcriptional heterogenetiy and does this correlate with CNV profiles

###########################################################
### Load sample-specific seurat objects
###########################################################

#perform unsupervised clustering based on gene expression data (optimized resolution)

fileName <- "MM1"

seurat.file <- paste0("./data/", fileName, "_seurat.rds")
seurat <- readRDS(seurat.file)

seurat <- StashIdent(seurat, save.name = "cluster.id")


###########################################################
### Load CNV subpopulation info
###########################################################

cell.groupings.file <- "./data/tsne_refined/icnv/icnv_cell_groupings.txt"

cell.groupings <- read.table(cell.groupings.file)
rownames(cell.groupings) <- cell.groupings$cell



###########################################################
### Subset CNV subpopulation info by sample and add to meta data
###########################################################

subset.cell.groupings <- cell.groupings[rownames(cell.groupings) %in% colnames(seurat), ]
subset.cell.groupings <- data.frame(subset.cell.groupings[, c("cell_group_name", "cell")])

#rename clones
og.names <- as.character(unique(subset.cell.groupings$cell_group_name))
new.names <- c(paste0("subpop_", 1:length(og.names)))

subset.cell.groupings$cell_group_name <- plyr::mapvalues(subset.cell.groupings$cell_group_name, 
                                                         from = og.names, 
                                                         to = new.names)


#add to seurat object

seurat <- subset(seurat, cells = rownames(subset.cell.groupings))

seurat <- AddMetaData(seurat, subset.cell.groupings)


###########################################################
### Create plots
###########################################################

meta <- seurat@meta.data





# --- UMAP - cluster colour

cent <- meta %>% group_by(cluster.id) %>% select(UMAP_1, UMAP_2) %>% summarize_all(median)

gg <- ggplot(meta, aes_string(x="UMAP_1", y="UMAP_2", col = "cluster.id"))
cluster <- gg + geom_point(size=1) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA, size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") +
  #scale_color_brewer(palette = "Set1") +
  labs(x="", y="", title = "Transcriptional Subpopulation") +
  geom_label_repel(aes(label = cluster.id), data = cent, show.legend = F, size = 3)  
cluster


# --- UMAP - CNV subpopulation colour

gg <- ggplot(meta, aes(UMAP_1, UMAP_2, color=cell_group_name))
cnv <- gg + geom_point(size=1) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA, size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  labs(x="", y="", title = "CNV Subpopulation") 
cnv


# --- Bar plot - xaxis is cluster, fill is CNV subpopulation colour 

bp <- ggplot(meta, aes(x=cluster.id, fill=cell_group_name)) +
  geom_bar(position = "fill") +
  #scale_fill_manual(values = c(colorPal_Neg[3], colorPal_Pos[3])) +
  labs(x="Cluster", y="Proportion") +
  theme(text = element_text(size = 14, family = "Helvetica", colour = "black"), 
        #axis.text.y = element_blank(), 
        #axis.ticks = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA, size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  labs(x="Transcriptional Subpopulation", y="Proportion") +
  scale_fill_brewer(palette = "Set1")
bp

# --- Plot together

cluster + 
cnv + 
bp

ggarrange(cluster, cnv, bp,
          ncol = 1)



# --- CNV-driven transcriptional subpopulations 

meta$cell_group_name <- as.character(meta$cell_group_name)
meta$cluster.id <- as.character(meta$cluster.id)

df <- data.frame(prop.table(x = table(meta$cell_group_name, meta$cluster.id),
                                       margin = 2))
df
```



