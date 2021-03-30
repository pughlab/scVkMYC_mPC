# Figure 2C-E

## Sample-wise DE anaylsis
#### Run DE analysis
```{r}

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_mPConly", "_seurat.rds")
dat <- readRDS(seurat.file)

Idents(dat) <- "Sample_ID"

DimPlot(dat, reduction = "tsne")


###########################################################
### Run DE analysis
###########################################################

samples <- c("ND1", "ND3", "ND4", "ND5",
             "AMG1", "AMG2", "AMG3",
             "MM1", "MM2", "MM3","MM4","MM5","MM6","MM7")


for (i in samples) {

  if (i %in% grep("ND", samples, value = T)) {
  
      dat_tmp1 <- subset(dat, subset = Cohort_ID == c("AMG", "MM"))
      dat_tmp2 <- subset(dat, subset = Sample_ID == i)
      dat_tmp <- subset(dat, cells = c(rownames(dat_tmp1@meta.data),
                                       rownames(dat_tmp2@meta.data)))
      
      cluster.markers.all <- FindMarkers(object = dat_tmp, 
                                      ident.1 = i, 
                                      min.pct = 0.3, 
                                      logfc.threshold = 0.5, 
                                      only.pos=F,
                                      pseudocount.use = 1/(nrow(dat_tmp@meta.data)))
      cluster.markers.all <- cluster.markers.all[cluster.markers.all$p_val_adj <= 0.05, ]
      
      cluster.markers.all$gene <- rownames(cluster.markers.all)
      cluster.markers.all <- cluster.markers.all[order(cluster.markers.all$avg_logFC, 
                                                       decreasing = T),]
      cluster.markers.all$Sample_ID <- i
      
      de.file <- paste0("./data/tsne_refined/SampleWise_DE_Analysis_Results_", i, ".csv")
      write.csv(cluster.markers.all, file = de.file)
  
  
  } 
  
  if (i %in% grep("AMG", samples, value = T)) {
    
       dat_tmp1 <- subset(dat, subset = Cohort_ID == c("ND", "MM"))
       dat_tmp2 <- subset(dat, subset = Sample_ID == i)
       dat_tmp <- subset(dat, cells = c(rownames(dat_tmp1@meta.data),
                                       rownames(dat_tmp2@meta.data)))
      
        cluster.markers.all <- FindMarkers(object = dat_tmp, 
                                        ident.1 = i, 
                                        min.pct = 0.3, 
                                        logfc.threshold = 0.5, 
                                        only.pos=F,
                                        pseudocount.use = 1/(nrow(dat_tmp@meta.data)))
        cluster.markers.all <- cluster.markers.all[cluster.markers.all$p_val_adj <= 0.05, ]
        
        cluster.markers.all$gene <- rownames(cluster.markers.all)
        cluster.markers.all <- cluster.markers.all[order(cluster.markers.all$avg_logFC, 
                                                         decreasing = T),]
        cluster.markers.all$Sample_ID <- i
        
        de.file <- paste0("./data/tsne_refined/SampleWise_DE_Analysis_Results_", i, ".csv")
        write.csv(cluster.markers.all, file = de.file) 
  
  } 
  
  if (i %in% grep("MM", samples, value = T)) {
    
       dat_tmp1 <- subset(dat, subset = Cohort_ID == c("ND", "AMG"))
       dat_tmp2 <- subset(dat, subset = Sample_ID == i)
       dat_tmp <- subset(dat, cells = c(rownames(dat_tmp1@meta.data),
                                       rownames(dat_tmp2@meta.data)))
      
        cluster.markers.all <- FindMarkers(object = dat_tmp, 
                                        ident.1 = i, 
                                        min.pct = 0.3, 
                                        logfc.threshold = 0.5, 
                                        only.pos=F,
                                        pseudocount.use = 1/(nrow(dat_tmp@meta.data)))
        cluster.markers.all <- cluster.markers.all[cluster.markers.all$p_val_adj <= 0.05, ]
        
        cluster.markers.all$gene <- rownames(cluster.markers.all)
        cluster.markers.all <- cluster.markers.all[order(cluster.markers.all$avg_logFC, 
                                                         decreasing = T),]
        cluster.markers.all$Sample_ID <- i
        
        de.file <- paste0("./data/tsne_refined/SampleWise_DE_Analysis_Results_", i, ".csv")
        write.csv(cluster.markers.all, file = de.file) 
  
  } 
    
}
  




```

#### Find overlaps 
##### ND
```{r}

###########################################################
### Load DE results 
###########################################################

# -- load de results

de.files <- grep("SampleWise_DE_Analysis_Results_ND", list.files("./data/tsne_refined"), value = T)

for (i in seq_along(de.files)) {
  
  if (i == as.integer(1)) {
    
    file <- paste0("./data/tsne_refined/", de.files[i])
    gs <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    
    }
  
    if (i > as.integer(1)){
      
    file <- paste0("./data/tsne_refined/", de.files[i])
    tmp <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    
    gs <- rbind(gs, tmp)
    
    }
}

dim(gs)

length(unique(gs$Sample_ID))

tmp.up <- gs[gs$avg_logFC > 0, ]
tmp.up <- tmp.up[tmp.up$pct.2 <= 0.3, ] #REMOVE PCT.2 > 0.3

tmp.down <- gs[gs$avg_logFC < 0, ]


table(tmp.up$Sample_ID)
table(tmp.down$Sample_ID)

# -- find overlaps 

write.csv(table(tmp.up$gene, tmp.up$Sample_ID),
          file = "./data/tsne_refined/SampleWise_DE_Analysis_ND_UP_table.csv")
write.csv(table(tmp.down$gene, tmp.down$Sample_ID),
          file = "./data/tsne_refined/SampleWise_DE_Analysis_ND_DOWN_table.csv")

#have to do this for formating reasons but don't need to change anything about the actiual csv file
df.up <- read.csv("./data/tsne_refined/SampleWise_DE_Analysis_ND_UP_table.csv", row.names = 1)
df.down <- read.csv("./data/tsne_refined/SampleWise_DE_Analysis_ND_DOWN_table.csv", row.names = 1)

#only genes in all groups 
overlaps_up <- df.up[rowSums( df.up != 0 ) >= length(unique(gs$Sample_ID)),] #;)
overlaps_down <- df.down[rowSums( df.down != 0 ) >= length(unique(gs$Sample_ID)),] #;)

#only use genes in 3 groups (USE THIS BECAUSE ND2 ONLY HAS 4 MPCs)
overlaps_up <- df.up[rowSums( df.up != 0 ) >= 3,] #;) #8
overlaps_down <- df.down[rowSums( df.down != 0 ) >= 3,] #;) #1

write.csv(overlaps_up, file = "./data/tsne_refined/SampleWise_DE_Analysis_ND_UP_genes.csv")
write.csv(overlaps_down, file = "./data/tsne_refined/SampleWise_DE_Analysis_ND_DOWN_genes.csv")

write.csv(tmp.up[tmp.up$gene %in% rownames(overlaps_up), ], file = "./data/tsne_refined/SampleWise_DE_Analysis_ND_UP_expression.csv")
write.csv(tmp.down[tmp.down$gene %in% rownames(overlaps_down), ], file = "./data/tsne_refined/SampleWise_DE_Analysis_ND_DOWN_expression.csv")




```


##### AMG
```{r}

###########################################################
### Load DE results 
###########################################################

# -- load de results

de.files <- grep("SampleWise_DE_Analysis_Results_AMG", list.files("./data/tsne_refined"), value = T)

for (i in seq_along(de.files)) {
  
  if (i == as.integer(1)) {
    
    file <- paste0("./data/tsne_refined/", de.files[i])
    gs <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    
    }
  
    if (i > as.integer(1)){
      
    file <- paste0("./data/tsne_refined/", de.files[i])
    tmp <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    
    gs <- rbind(gs, tmp)
    
    }
}

dim(gs)

length(unique(gs$Sample_ID))

tmp.up <- gs[gs$avg_logFC > 0, ]
tmp.up <- tmp.up[tmp.up$pct.2 <= 0.3, ] #REMOVE PCT.2 > 0.4

tmp.down <- gs[gs$avg_logFC < 0, ]

table(tmp.up$Sample_ID)
table(tmp.down$Sample_ID)

# -- find overlaps 

write.csv(table(tmp.up$gene, tmp.up$Sample_ID),
          file = "./data/tsne_refined/SampleWise_DE_Analysis_AMG_UP_table.csv")
write.csv(table(tmp.down$gene, tmp.down$Sample_ID),
          file = "./data/tsne_refined/SampleWise_DE_Analysis_AMG_DOWN_table.csv")

#have to do this for formating reasons but don't need to change anything about the actiual csv file
df.up <- read.csv("./data/tsne_refined/SampleWise_DE_Analysis_AMG_UP_table.csv", row.names = 1)
df.down <- read.csv("./data/tsne_refined/SampleWise_DE_Analysis_AMG_DOWN_table.csv", row.names = 1)

#only genes in all groups 
overlaps_up <- df.up[rowSums( df.up != 0 ) >= length(unique(gs$Sample_ID)),] #;)
overlaps_down <- df.down[rowSums( df.down != 0 ) >= length(unique(gs$Sample_ID)),] #;)
# 
# #only use genes in 3 groups (USE THIS BECAUSE ND2 ONLY HAS 4 MPCs)
# overlaps_up <- df.up[rowSums( df.up != 0 ) >= 3,] #;) #8
# overlaps_down <- df.down[rowSums( df.down != 0 ) >= 3,] #;) #1

write.csv(overlaps_up, file = "./data/tsne_refined/SampleWise_DE_Analysis_AMG_UP_genes.csv")
write.csv(overlaps_down, file = "./data/tsne_refined/SampleWise_DE_Analysis_AMG_DOWN_genes.csv")

write.csv(tmp.up[tmp.up$gene %in% rownames(overlaps_up), ], file = "./data/tsne_refined/SampleWise_DE_Analysis_AMG_UP_expression.csv")
write.csv(tmp.down[tmp.down$gene %in% rownames(overlaps_down), ], file = "./data/tsne_refined/SampleWise_DE_Analysis_AMG_DOWN_expression.csv")




```


##### MM
```{r}

###########################################################
### Load DE results 
###########################################################

# -- load de results

de.files <- grep("SampleWise_DE_Analysis_Results_MM", list.files("./data/tsne_refined"), value = T)

for (i in seq_along(de.files)) {
  
  if (i == as.integer(1)) {
    
    file <- paste0("./data/tsne_refined/", de.files[i])
    gs <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    
    }
  
    if (i > as.integer(1)){
      
    file <- paste0("./data/tsne_refined/", de.files[i])
    tmp <- read.csv(file, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
    
    gs <- rbind(gs, tmp)
    
    }
}

dim(gs)

length(unique(gs$Sample_ID))

tmp.up <- gs[gs$avg_logFC > 0, ]
tmp.up <- tmp.up[tmp.up$pct.2 <= 0.3, ] #REMOVE PCT.2 > 0.4

tmp.down <- gs[gs$avg_logFC < 0, ]

table(tmp.up$Sample_ID)
table(tmp.down$Sample_ID)

# -- find overlaps 

write.csv(table(tmp.up$gene, tmp.up$Sample_ID),
          file = "./data/tsne_refined/SampleWise_DE_Analysis_MM_UP_table.csv")
write.csv(table(tmp.down$gene, tmp.down$Sample_ID),
          file = "./data/tsne_refined/SampleWise_DE_Analysis_MM_DOWN_table.csv")

#have to do this for formating reasons but don't need to change anything about the actiual csv file
df.up <- read.csv("./data/tsne_refined/SampleWise_DE_Analysis_MM_UP_table.csv", row.names = 1)
df.down <- read.csv("./data/tsne_refined/SampleWise_DE_Analysis_MM_DOWN_table.csv", row.names = 1)

#only genes in all groups 
# overlaps_up <- df.up[rowSums( df.up != 0 ) >= length(unique(gs$Sample_ID)),] #;)
# overlaps_down <- df.down[rowSums( df.down != 0 ) >= length(unique(gs$Sample_ID)),] #;)

#USE THIS
overlaps_up <- df.up[rowSums( df.up != 0 ) >= 7,] #;) #no overlapping in 5,6,7 (suggests diversification of tumours - inter-tumoural)
overlaps_down <- df.down[rowSums( df.down != 0 ) >= 7,] #;) #8

write.csv(overlaps_up, file = "./data/tsne_refined/SampleWise_DE_Analysis_MM_UP_genes.csv")
write.csv(overlaps_down, file = "./data/tsne_refined/SampleWise_DE_Analysis_MM_DOWN_genes.csv")

write.csv(tmp.up[tmp.up$gene %in% rownames(overlaps_up), ], file = "./data/tsne_refined/SampleWise_DE_Analysis_MM_UP_expression.csv")
write.csv(tmp.down[tmp.down$gene %in% rownames(overlaps_down), ], file = "./data/tsne_refined/SampleWise_DE_Analysis_MM_DOWN_expression.csv")


```


### Barplots for sample-wise DE - with Cont PCs only

```{r}

###########################################################
### Load seurat object 
###########################################################

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_PConly", "_seurat.rds")
seurat <- readRDS(seurat.file)
npc <- WhichCells(seurat, idents = "nPC")
npc <- grep("Cont", npc, value = T)

mpc <-  WhichCells(seurat, idents = "mPC")
seurat <- StashIdent(seurat, save.name = "status")

seurat@meta.data$Sample_ID_2 <- paste(seurat@meta.data$Sample_ID, 
                                      seurat@meta.data$status, sep = "_")


Idents(seurat) <- "Sample_ID_2"
DimPlot(seurat)

seurat <- subset(seurat, cells = c(npc, mpc))
seurat <- subset(seurat, subset = Sample_ID %in% c("ND2", "ND3"), invert = T)
DimPlot(seurat)



###########################################################
### Fetch Genes for each plot 
###########################################################

# Differentially expressed genes

dir <- "./data/tsne_refined/SampleWise_DE_Analysis_"
goi.full <- c(rownames(read.csv(paste0(dir, "ND_UP_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "ND_DOWN_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "AMG_UP_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "AMG_DOWN_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "MM_UP_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "MM_DOWN_genes.csv"), row.names = 1))
         )

goi <- unique(goi.full)



###########################################################
### Compute average expression of each gene for each sample 
###########################################################

# Re-scale data
seurat <- ScaleData(seurat, vars.to.regress = "percent.mt", 
                    features = goi)

avg.exp.sample <- data.frame(AverageExpression(seurat,
                                               slot = "scale.data",
                                               features = goi))
colnames(avg.exp.sample) <- gsub("RNA.", "", colnames(avg.exp.sample))

module_scores <- data.frame(t(avg.exp.sample))

#Sort by sample_ID
module_scores$Sample_ID <- rownames(module_scores)


#Add cohort info 

module_scores$Cohort_ID <- rownames(module_scores)

module_scores[rownames(module_scores) %in% 
                grep("AMG",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "int-MM"

module_scores[rownames(module_scores) %in% 
                grep("^MM",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "active-MM"

module_scores[rownames(module_scores) %in% 
                grep("ND",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "early-MM"

module_scores[rownames(module_scores) %in% 
                grep("nPC",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "nPC"



module_scores$Cohort_ID <- factor(module_scores$Cohort_ID, 
                                  levels = c("nPC", "early-MM", "int-MM", "active-MM"))

###########################################################
### Line barplot with overlapid line plot including stats
###########################################################

melted <- melt(module_scores)

colnames(melted) <- gsub("variable", "gene", colnames(melted))

# Set up stats (https://www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets/)

library(ggpubr)
library(rstatix)
library(RColorBrewer)
stat.test <- melted %>%
  group_by(gene) %>%
  t_test(value ~ Cohort_ID, p.adjust.method = "bonferroni") #t_test or wilcox_test
stat.test

stat.test <- stat.test %>% add_y_position()



ggplot(melted, aes(x=Cohort_ID, y=value)) +
  geom_boxplot(aes(fill = Cohort_ID), outlier.shape = NA) +
  facet_wrap(~gene,   scales = "free") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_base() +
  scale_fill_manual(values = brewer.pal(4, "Set2")[1:4])






```

## Line plots
```{r}

## THIS WAS USED IN MANUSCRIPT 

###########################################################
### Load seurat object 
###########################################################

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_PConly", "_seurat.rds")
seurat <- readRDS(seurat.file)
npc <- WhichCells(seurat, idents = "nPC")
npc <- grep("Cont", npc, value = T)

mpc <-  WhichCells(seurat, idents = "mPC")
seurat <- StashIdent(seurat, save.name = "status")

seurat@meta.data$Sample_ID_2 <- paste(seurat@meta.data$Sample_ID, 
                                      seurat@meta.data$status, sep = "_")


Idents(seurat) <- "Sample_ID_2"
DimPlot(seurat)

seurat <- subset(seurat, cells = c(npc, mpc))
seurat <- subset(seurat, subset = Sample_ID %in% c("ND2", "ND3"), invert = T)
DimPlot(seurat)



###########################################################
### Fetch Genes for each plot 
###########################################################

# Differentially expressed genes

dir <- "./data/tsne_refined/SampleWise_DE_Analysis_"
goi.full <- c(rownames(read.csv(paste0(dir, "ND_UP_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "ND_DOWN_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "AMG_UP_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "AMG_DOWN_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "MM_UP_genes.csv"), row.names = 1)),
         rownames(read.csv(paste0(dir, "MM_DOWN_genes.csv"), row.names = 1))
         )

goi <- unique(goi.full)



###########################################################
### Compute average expression of each gene for each sample 
###########################################################

# Re-scale data
seurat <- ScaleData(seurat, vars.to.regress = "percent.mt", 
                    features = goi)

avg.exp.sample <- data.frame(AverageExpression(seurat,
                                               slot = "scale.data",
                                               features = goi))
colnames(avg.exp.sample) <- gsub("RNA.", "", colnames(avg.exp.sample))

module_scores <- data.frame(t(avg.exp.sample))

#Sort by sample_ID
module_scores$Sample_ID <- rownames(module_scores)


#Add cohort info 

module_scores$Cohort_ID <- rownames(module_scores)

module_scores[rownames(module_scores) %in% 
                grep("AMG",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "int-MM"

module_scores[rownames(module_scores) %in% 
                grep("MM",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "active-MM"

module_scores[rownames(module_scores) %in% 
                grep("ND",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "early-MM"

module_scores[rownames(module_scores) %in% 
                grep("nPC",
                     rownames(module_scores), 
                     value = T), 
              "Cohort_ID"] <- "nPC"



module_scores$Cohort_ID <- factor(module_scores$Cohort_ID, 
                                  levels = c("nPC", "early-MM", "int-MM", "active-MM"))




###########################################################
### Prep data for plotting
###########################################################

module_scores$Sample_ID <- NULL

library(Rmisc)
tgc <- summarySE(melted, measurevar="value", groupvars=c("Cohort_ID","gene"))

#sig.genes (based on calculations from "Barplots for sample-wise DE - with Cont PCs only")
sig.genes <- as.character(unique(tgc$gene)[!unique(tgc$gene) %in% 
                                    c("Cotl1", "Cxcl2", "G0s2", "Gadd45a", "Gm10800", "Gm26870",
                                      "Gm6563", "H2.Aa", "Laptm5", "Mgst1", "Ms4a3", "Rgs2", "Trem3", 
                                      "H2.Ab1", "H2.Eb1", "Lars2", "mt.Atp8")])

tgc <- tgc[tgc$gene %in% sig.genes, ]
tgc$gene <- factor(as.character(tgc$gene),
                       levels = c(c("Gm2a", "Myl4", 
                                    "Il5ra", "Mgmt", "Tsc22d1",
                                    sig.genes[!sig.genes %in%  c("Gm2a", "Myl4", 
                                    "Il5ra", "Mgmt", "Tsc22d1")])
                                  )
                       )


# ggplot(tgc, aes(x=Cohort_ID, y=value, group=gene)) + 
#     geom_errorbar(aes(ymin=value-se, ymax=value+se, color=gene), width=.1) +
#     geom_line(aes(color=gene)) +
#   facet_wrap(~gene, ncol = 6, scales = "free_y") +
#   theme_calc() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.position = "none")


# for sig genes (Fig 2c)

ggplot(tgc, aes(x=Cohort_ID, y=value, group=gene)) + 
    geom_errorbar(aes(ymin=value-se, ymax=value+se, color=Cohort_ID), width=.1) +
    geom_line(color = "grey", size=0.5,) +
  geom_point(aes(color=Cohort_ID)) +
  facet_wrap(~gene, ncol = 6, scales = "free_y") +
  theme_base() +
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        panel.border = element_rect(color="black"), 
        legend.position = "bottom", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Set2")






```


