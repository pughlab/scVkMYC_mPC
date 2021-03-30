# Figure 2A

## DEG between cohorts
```{r}

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_PConly", "_seurat.rds")

seurat <- readRDS(seurat.file)

mm <- subset(seurat, subset = Cohort_ID == "MM") 
mm_mP <- subset(mm, idents = c("nPC"), invert = TRUE)

amg <- subset(seurat, subset = Cohort_ID == "AMG") 
amg_mP <- subset(amg, idents = c("nPC"), invert = TRUE)

nd <- subset(seurat, subset = Cohort_ID == "ND") 
nd_mP <- subset(nd, idents = c("nPC"), invert = TRUE)

Idents(seurat, cells = rownames(mm_mP@meta.data)) <- "MM_mPC"
Idents(seurat, cells = rownames(amg_mP@meta.data)) <- "AMG_mPC"
Idents(seurat, cells = rownames(nd_mP@meta.data)) <- "ND_mPC"

DimPlot(seurat, reduction = "tsne")

## ----- cluster.markers.AMG_mPC


## - AMG mPC vs. all other PCs
cluster.markers.AMG_mPC <- FindMarkers(object = seurat, 
                                       ident.1 = "AMG_mPC", 
                                       ident.2 = c("nPC", "MM_mPC", "ND_mPC"),
                                       test.use = "wilcox",
                                       min.pct = 0.3, 
                                       logfc.threshold = 0.5, 
                                       only.pos=F,
                                       pseudocount.use = 1/nrow(seurat@meta.data))


cluster.markers.AMG_mPC <- cluster.markers.AMG_mPC[cluster.markers.AMG_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.AMG_mPC, file = "./data/tsne_refined/DEanalysis_AMGmPC_vs_allPC.csv")

## - AMG mPC vs. nPCs
cluster.markers.AMG_mPC <- FindMarkers(object = seurat, 
                                       ident.1 = "AMG_mPC", 
                                       ident.2 = c("nPC"),
                                       test.use = "wilcox",
                                       min.pct = 0.3, 
                                       logfc.threshold = 0.5, 
                                       only.pos=F,
                                       pseudocount.use = 1/nrow(seurat@meta.data))


cluster.markers.AMG_mPC <- cluster.markers.AMG_mPC[cluster.markers.AMG_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.AMG_mPC, file = "./data/tsne_refined/DEanalysis_AMGmPC_vs_nPC.csv")


## - AMG mPC vs. all other mPCs
cluster.markers.AMG_mPC <- FindMarkers(object = seurat, 
                                       ident.1 = "AMG_mPC", 
                                       ident.2 = c("MM_mPC", "ND_mPC"),
                                       test.use = "wilcox",
                                       min.pct = 0.3, 
                                       logfc.threshold = 0.5, 
                                       only.pos=F,
                                       pseudocount.use = 1/nrow(seurat@meta.data))


cluster.markers.AMG_mPC <- cluster.markers.AMG_mPC[cluster.markers.AMG_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.AMG_mPC, file = "./data/tsne_refined/DEanalysis_AMGmPC_vs_mPC.csv")



## ----- cluster.markers.ND_mPC

## - ND mPC vs. all other PCs
cluster.markers.ND_mPC <- FindMarkers(object = seurat, 
                                      ident.1 = "ND_mPC", 
                                      ident.2 = c("nPC","MM_mPC", "AMG_mPC"),
                                      min.pct = 0.3, 
                                      logfc.threshold = 0.5, 
                                      only.pos=F,
                                      pseudocount.use = 1/nrow(seurat@meta.data))

cluster.markers.ND_mPC <- cluster.markers.ND_mPC[cluster.markers.ND_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.ND_mPC, file = "./data/tsne_refined/DEanalysis_NDmPC_vs_allPC.csv")

## - ND mPC vs. nPCs
cluster.markers.ND_mPC <- FindMarkers(object = seurat, 
                                       ident.1 = "ND_mPC", 
                                       ident.2 = c("nPC"),
                                       test.use = "wilcox",
                                       min.pct = 0.3, 
                                       logfc.threshold = 0.5, 
                                       only.pos=F,
                                       pseudocount.use = 1/nrow(seurat@meta.data))

cluster.markers.ND_mPC <- cluster.markers.ND_mPC[cluster.markers.ND_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.ND_mPC, file = "./data/tsne_refined/DEanalysis_NDmPC_vs_nPC.csv")




## - ND mPC vs. all other mPCs
cluster.markers.ND_mPC <- FindMarkers(object = seurat, 
                                       ident.1 = "ND_mPC", 
                                       ident.2 = c("MM_mPC", "AMG_mPC"),
                                       test.use = "wilcox",
                                       min.pct = 0.3, 
                                       logfc.threshold = 0.5, 
                                       only.pos=F,
                                       pseudocount.use = 1/nrow(seurat@meta.data))


cluster.markers.ND_mPC <- cluster.markers.ND_mPC[cluster.markers.ND_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.ND_mPC, file = "./data/tsne_refined/DEanalysis_NDmPC_vs_mPC.csv")




## ----- cluster.markers.MM_mPC

## - MM mPC vs. all other PCs
cluster.markers.MM_mPC <- FindMarkers(object = seurat, 
                                      ident.1 = "MM_mPC", 
                                      ident.2 = c("nPC","ND_mPC", "AMG_mPC"),
                                      min.pct = 0.3, 
                                      logfc.threshold = 0.5, 
                                      only.pos=F,
                                      pseudocount.use = 1/nrow(seurat@meta.data))


cluster.markers.MM_mPC <- cluster.markers.MM_mPC[cluster.markers.MM_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.MM_mPC, file = "./data/tsne_refined/DEanalysis_MMmPC_vs_allPC.csv")

## - ND mPC vs. all other PCs
cluster.markers.MM_mPC <- FindMarkers(object = seurat, 
                                      ident.1 = "MM_mPC", 
                                      ident.2 = "nPC",
                                      min.pct = 0.3, 
                                      logfc.threshold = 0.5, 
                                      only.pos=F,
                                      pseudocount.use = 1/nrow(seurat@meta.data))

cluster.markers.MM_mPC <- cluster.markers.MM_mPC[cluster.markers.MM_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.MM_mPC, file = "./data/tsne_refined/DEanalysis_MMmPC_vs_nPC.csv")


## - ND mPC vs. all other mPCs
cluster.markers.MM_mPC <- FindMarkers(object = seurat, 
                                       ident.1 = "MM_mPC", 
                                       ident.2 = c("ND_mPC", "AMG_mPC"),
                                       test.use = "wilcox",
                                       min.pct = 0.3, 
                                       logfc.threshold = 0.5, 
                                       only.pos=F,
                                       pseudocount.use = 1/nrow(seurat@meta.data))


cluster.markers.MM_mPC <- cluster.markers.MM_mPC[cluster.markers.MM_mPC$p_val_adj <= 0.05, ]
write.csv(cluster.markers.MM_mPC, file = "./data/tsne_refined/DEanalysis_MMmPC_vs_mPC.csv")



```

## Overlapping genes (core) 

```{r}

#load DE results for mPC of disease group vs. nPCs 

disease <- c("ND", "AMG", "MM")

for (i in disease) { 
  
  if (i == "ND") {
    
    de.file <- paste0("./data/tsne_refined/DEanalysis_", i, "mPC_vs_nPC.csv")
    de.results <- read.csv(de.file, row.names = 1) 
    de.results$Cohort_ID <- i
    de.results$gene <- rownames(de.results)
    
    
    de.results.up <- de.results[de.results$avg_logFC > 0, ]
    de.results.down <- de.results[de.results$avg_logFC < 0, ]
    
  } else {
    
    de.file <- paste0("./data/tsne_refined/DEanalysis_", i, "mPC_vs_nPC.csv")
    tmp <- read.csv(de.file, row.names = 1) 
    tmp$Cohort_ID <- i
    tmp$gene <- rownames(tmp)
    
    tmp.up <- tmp[tmp$avg_logFC > 0, ]
    tmp.down <- tmp[tmp$avg_logFC < 0, ]
    
    de.results.up <- rbind(de.results.up, tmp.up)
    de.results.down <- rbind(de.results.down, tmp.down)

    
  }
}


# Find overlaps of UP
write.csv(table(de.results.up$gene, de.results.up$Cohort_ID),
          file = "./data/tsne_refined/deresultsCORE_UP_table.csv")

#have to do this for formating reasons but don't need to change anything about the actiual csv file
df <- read.csv("./data/tsne_refined/deresultsCORE_UP_table.csv", row.names = 1)

#only genes in all (3) groups 
overlaps <- df[rowSums( df != 0 ) >= 3,] #;)

write.csv(overlaps, file = "./data/tsne_refined/deresultsCORE_UP_Genes.csv")

write.csv(de.results.up[de.results.up$gene %in% rownames(overlaps), ], file = "./data/tsne_refined/deresultsCORE_UP_Expression.csv")


# Find overlaps of DOWN

write.csv(table(de.results.down$gene, de.results.down$Cohort_ID),
          file = "./data/tsne_refined/deresultsCORE_DOWN_table.csv")

#have to do this for formating reasons but don't need to change anything about the actiual csv file
df <- read.csv("./data/tsne_refined/deresultsCORE_DOWN_table.csv", row.names = 1)

#only genes in all (3) groups 
overlaps <- df[rowSums( df != 0 ) >= 3,] #;)
write.csv(overlaps, file = "./data/tsne_refined/deresultsCORE_DOWN_Genes.csv")

write.csv(de.results.down[de.results.down$gene %in% rownames(overlaps), ], file = "./data/tsne_refined/deresultsCORE_DOWN_Expression.csv")

```

### Heatmap

```{r}

###########################################################
### Load seurat object 
###########################################################

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_PConly", "_seurat.rds")
seurat <- readRDS(seurat.file)
seurat@meta.data$malig_status <- Idents(seurat)


###########################################################
### Fetch Genes for each plot 
###########################################################

# Differentially expressed genes

core.up <- read.csv("./data/tsne_refined/deresultsCORE_UP_Genes.csv", row.names = 1)
core.down <- read.csv("./data/tsne_refined/deresultsCORE_DOWN_Genes.csv", row.names = 1)

goi <- c(rownames(core.up), rownames(core.down))

# Re-scale data

seurat <- ScaleData(seurat, vars.to.regress = "percent.mt", 
                    features = goi)

# ---------------Downsample

#select 100 nPCs
npc <- subset(seurat, idents = "nPC") #576
npc_ds <- subset(npc, downsample = 100)

#select mPC - select 100 from each VkMYC cohort
mpc <- subset(seurat, idents = "mPC")
Idents(mpc) <- "Cohort_ID"
mpc_ds <- subset(mpc, downsample = 100)


seurat_ds <- subset(seurat, cells = c(rownames(npc_ds@meta.data), 
                                      rownames(mpc_ds@meta.data)))



genes_scaled <- FetchData(seurat_ds, vars = goi, slot = "scale.data")

#set min and max

genes_scaled <- MinMax(data = genes_scaled, min = -2, max = 2)


# fix sample and cohort samples

library(magrittr)

seurat_ds@meta.data$Cohort_ID %<>%
  gsub("Cont", "Control", .) %>%
   gsub("MM", "active-MM", .) %>%
  gsub("ND", "early-MM", .) %>%
  gsub("AMG", "int-MM", .)
  

seurat_ds@meta.data$Sample_ID %<>%
  gsub("ND", "EMM", .) %>%
  gsub("AMG", "IMM", .) %>%
  gsub("^MM", "AMM", .) 

seurat_ds@meta.data$Sample_ID <- factor(as.character(seurat_ds@meta.data$Sample_ID),
                                        levels = c("Cont1", "Cont2", "Cont3",
                                                   "EMM1", "EMM2", "EMM3", "EMM4", "EMM5",
                                                   "IMM1", "IMM2", "IMM3",
                                                   "AMM1", "AMM2", "AMM3","AMM4","AMM5","AMM6","AMM7"))


seurat_ds@meta.data$Cohort_ID <- factor(as.character(seurat_ds@meta.data$Cohort_ID),
                                        levels = c("Control", "early-MM",
                                                   "int-MM", "active-MM"))


#colours 

# heatmap_colours <- rev(brewer.pal(n = 7, name = "RdYlBu"))
 

gg_color_hue <- function(n) { #Nina's function
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pc_col <- c("grey", "#A75B5D")
cohort_col <- brewer.pal(4, "Set2")
sample_col <- gg_color_hue(18)

  
anno_colour <- list(PCstatus = c("nPC" = pc_col[1],
                                 "mPC" = pc_col[2]),
                    Cohort_ID = c("Control" = cohort_col[1],
                                  "early-MM" = cohort_col[2],
                                  "int-MM" = cohort_col[3],
                                  "active-MM" = cohort_col[4]),
                    Sample_ID = c("Cont1" = sample_col[1],
                                  "Cont2" = sample_col[2],
                                  "Cont3" = sample_col[3],
                                  "EMM1" = sample_col[4],
                                  "EMM3" = sample_col[6],
                                  "EMM4" = sample_col[7],
                                  "EMM5" = sample_col[8],
                                  "IMM1" = sample_col[9],
                                  "IMM2" = sample_col[10],
                                  "IMM3" = sample_col[11],
                                  "AMM1" = sample_col[12],
                                  "AMM2" = sample_col[13],
                                  "AMM3" = sample_col[14],
                                  "AMM4" = sample_col[15],
                                  "AMM5" = sample_col[16],
                                  "AMM6" = sample_col[17],
                                  "AMM7" = sample_col[18])
                                        )




#annotation bar 

annoFull <- seurat_ds@meta.data


annoFull <- annoFull %>% arrange(Cohort_ID) %>% arrange(malig_status) 

head(annoFull)

library(ComplexHeatmap)
column_anno <- HeatmapAnnotation(PCstatus = annoFull$malig_status,
                                 Cohort_ID = annoFull$Cohort_ID,
                                 Sample_ID = annoFull$Sample_ID,
                                 #gp = gpar(col = "black"),
                                 which = "column",
                                 border = TRUE,
                                 col = anno_colour
                                 )



#need to order genes_scaled
genes_scaled$Row.names <- rownames(genes_scaled)

genes_scaled$Row.names <- factor(genes_scaled$Row.names,
                                 levels = rownames(annoFull))
genes_scaled <- genes_scaled %>% arrange(Row.names)
genes_scaled$Row.names <- NULL




#plot

genes_scaled <- genes_scaled[, !colnames(genes_scaled) == "MYC"]

Heatmap(t(genes_scaled), name = "Scaled", 
        top_annotation = column_anno,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        border = TRUE, 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        #col = heatmap_colours,
        #column_order = 200:1,
        column_split = c(rep("active-MM", 100), rep("int-MM", 100), 
                         rep("early-MM", 100), rep("nPC", 100)),
        cluster_column_slices = FALSE,
        row_split = c(rep("A", nrow(core.up)-1),
                      rep("B", nrow(core.down)))
        #width = unit(8, "cm")
        )






```

# Figure 2B

### GSEA results 

```{r}

###########################################################
### GSEA
###########################################################

#using msigdb, C2, C6, H
#top20, FDR 0.05, mouse

#save - reformat and save as csv 
#reload data back into R


###########################################################
### Plotting GSEA results
###########################################################

df <- read.csv("./data/tsne_refined/deresultsCORE_GSEA.csv")


library(stringr)

#make enrichment for core terms negative

gg <- ggplot(df, aes(y=k.K.Direction,x=reorder(Gene.Set.Name,k.K.Direction), fill=-log10(FDR.q.value)))
gg + geom_bar(stat="identity", colour = "black") + 
  coord_flip() + 
  labs(y="Reactions Ratio", x="Reactome Term", fill = "-log10 FDR") +
  scale_fill_viridis_c(option="plasma",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() + theme_base() +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dashed", size = 0.25),
        panel.grid.minor.x = element_line(colour = "grey", linetype = "dashed", size = 0.25),
        axis.text = element_text(size = 13), legend.position = c(-0.9, 0.55),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust=0),
        panel.border = element_rect(linetype = "solid",
                                                fill = NA, size = 0.75),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1.5,"cm") 
        ) 
  #scale_x_discrete(position = "top") 
  #coord_flip()







```

