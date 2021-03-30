# Fig 3 Supplemental
## Reference heatmap
```{r}

#after running iCNV 

setwd(outputDir)
icnv.dir <- "./data/tsne_refined/icnv/scVKMYC_TSNErefined_PConly_HMMM/"

observations <- read.csv(paste0(icnv.dir, "infercnv.references.txt"), sep = " ")

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
### Set up colours 
###########################################################

#Metadata annotations 
gg_color_hue <- function(n) { #Nina's function
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# cohort.col <- brewer.pal(4, "Set2")[2:4]
# sample.col <- gg_color_hue(18)[c(4,6:18)]


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
### Set up meta data
###########################################################

# seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_mPConly", "_seurat.rds")
# dat <- readRDS(seurat.file)
# 
# # fix sample and cohort samples
# 
# library(magrittr)
# 
# dat@meta.data$Cohort_ID %<>%
#    gsub("MM", "active-MM", .) %>%
#   gsub("ND", "early-MM", .) %>%
#   gsub("AMG", "int-MM", .)
#   
# 
# dat@meta.data$Sample_ID %<>%
#   gsub("ND", "EMM", .) %>%
#   gsub("AMG", "IMM", .) %>%
#   gsub("^MM", "AMM", .) 
# 
# dat@meta.data$Sample_ID <- factor(as.character(dat@meta.data$Sample_ID),
#                                         levels = c(
#                                                    "EMM1", "EMM3", "EMM4", "EMM5",
#                                                    "IMM1", "IMM2", "IMM3",
#                                                    "AMM1", "AMM2", "AMM3","AMM4","AMM5","AMM6","AMM7"))
# 
# 
# dat@meta.data$Cohort_ID <- factor(as.character(dat@meta.data$Cohort_ID),
#                                         levels = c("early-MM",
#                                                    "int-MM", "active-MM"))
# 
# 
# meta <- dat@meta.data[rownames(dat@meta.data) %in% rownames(observations_mat), ]
# 
# #re-order to match observations_mat
# meta$Row.names <- rownames(meta)
# meta$Row.names <- factor(meta$Row.names,
#                          levels = rownames(observations_mat))
# meta <- meta %>% arrange(Row.names)
# meta$Row.names <- NULL
# 


###########################################################
### Specify annotation colours 
###########################################################

anno_colours_column <- list(Chromosome = chr.colours)
names(anno_colours_column$Chromosome) <- levels(ordered_genes$V2)

# anno_colours_rows <- list(Cohort_ID = cohort.col,
#                           Sample_ID = sample.col)
# 
# names(anno_colours_rows$Cohort_ID) <- levels(meta$Cohort_ID)
# names(anno_colours_rows$Sample_ID) <- levels(meta$Sample_ID)
# 



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

# row_anotations <- HeatmapAnnotation(Sample_ID = meta$Sample_ID,
#                                     Cohort_ID = meta$Cohort_ID, 
#                                     which = "row", 
#                                     border = TRUE,
#                                     col = anno_colours_rows, 
#                                     annotation_legend_param = list(Sample_ID = list(direction = "horizontal"),
#                                                                    Cohort_ID = list(direction = "horizontal")))


###########################################################
### Plot heatmap, 
###########################################################


Heatmap(observations_mat, 
        name = "matrix_9",
                top_annotation = column_annotations,
                #right_annotation = row_anotations,
                #show_row_dend = TRUE,
                #row_split = 8,
                show_column_dend = FALSE,
                border = TRUE, 
                row_names_gp = gpar(fontsize = 12),
                cluster_rows = FALSE,
                #clustering_distance_rows = "euclidean",
                #clustering_method_rows = "ward.D",
                cluster_columns = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                col = heatmap.cols,
                cluster_column_slices = FALSE,
                column_split = c(ordered_genes$V2),
                row_dend_width = unit(2, "cm")
                )







```

# Chr5 reactome
```{r}

cell.groupings.file <- "./data/tsne_refined/icnv/icnv_cell_groupings.txt"

cell.groupings <- read.table(cell.groupings.file)
head(cell.groupings)
#cell.groupings <- cell.groupings[cell.groupings$sample == "MM1", ]



cell.groupings.file <- "./data/tsne_refined/icnv/icnv_pred_cnv_groupings.txt"

pred_cnv <- read.table(cell.groupings.file)

#pred_cnv <- pred_cnv[pred_cnv$sample == "MM1", ]
pred_cnv$chr_state <- paste(pred_cnv$chr, pred_cnv$state, sep = "_")
head(pred_cnv)


#select cell_group_names with chr5 state 1

chr5del_group <- pred_cnv[pred_cnv$chr_state == "chr5_1", ]
chr5del_group <- as.character(unique(chr5del_group$cell_group_name))

cell.groupings[cell.groupings$cell_group_name %in% chr5del_group, "chr5status"] <- "del"
cell.groupings[is.na(cell.groupings$chr5status), "chr5status"] <- "wt"

head(cell.groupings)




##Chr5 --> Reactome

# --- Load seurat data

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_mPConly", "_seurat.rds")
subset <- readRDS(seurat.file)

DimPlot(subset)


# --- Add chr5 status to metadata


cell.groupings <- cell.groupings[rownames(cell.groupings) %in% rownames(subset@meta.data), ]
subset <- AddMetaData(subset, cell.groupings)

# --- Set ident to chr5 status and run DGE analysis


Idents(subset) <- "chr5status"

DimPlot(subset)

cluster.markers.all <- FindMarkers(object = subset, 
                                   ident.1 = "del", 
                                   ident.2 = "wt",
                                        min.pct = 0.3, 
                                        logfc.threshold = 1.0, 
                                        only.pos=F,
                                        pseudocount.use = 1/(nrow(subset@meta.data)))

cluster.markers.all <- cluster.markers.all[cluster.markers.all$p_val_adj <= 0.05, ]

# --- Save DGE analysis and use for reactome enrichment analysis 

write.csv(cluster.markers.all, file = "./data/scVKMYC_TSNErefined_DEG_chr5vsWT_mPCs.csv")


#bring it back in
#library(Hmisc)
library(dplyr)
library(reshape2)
library(stringr)

reactome.files <- "chr5_up_reactome"

file <- paste0("./data/tsnerefined_chr5_up_reactome.csv")
gs <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")

dim(gs)
#[1] 206  15


#remove entries with FDR > 0.05
gs <- gs[gs$Entities.FDR <= 0.05, ] 

dim(gs)


#remove entries with only 1 gene (Submitted.entities.found)
gs <- gs[grep(";", gs$Submitted.entities.found), ] 

dim(gs)

library(stringr)
tmp <- data.frame(gs$Submitted.entities.found)
gs$count <- apply(tmp, 1, function(x) str_count(x, pattern = ";")) + 1
rm(tmp)

dim(gs)
#[1] 70 15




# Plot reactions ratio with fdr
 
library(stringr)

gs_20 <- gs %>% arrange(-Reactions.ratio) %>% head(n= 20)
gs_20$newx = str_wrap(gs_20$Pathway.name, width = 50)


gg <- ggplot(gs_20, aes(y=Reactions.ratio,x=reorder(newx,Reactions.ratio), fill=-log10(Entities.FDR)))
gg + geom_bar(stat="identity", colour = "black") + 
  coord_flip() + 
  labs(y="Enrichment Score", x="Reactome Term", fill = "-log10 FDR") +
  scale_fill_viridis_c(option="plasma",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +   theme_base() +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", size = 0.25),
        panel.grid.minor.x = element_line(colour = "grey", linetype = "dashed", size = 0.25),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 20), 
        legend.position = c(0.9, 0.25),
        axis.title = element_text(size = 15),
        panel.border = element_rect(linetype = "solid",
                                                fill = NA, size = 0.75),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1.5,"cm") 
        ) 
  #theme(legend.position = "none") 
        


```

