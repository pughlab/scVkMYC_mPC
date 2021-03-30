
### Figure 1B. B/PC TSNE

```{r}

###########################################################
### Load input data
###########################################################


seuratObj <- "/Volumes/Samsung_T5/OneDrive - UHN/Documents/scVKMYC/analysis/cohort/lineages/scVKMYC_cohort_lineages/B_PC/scVKMYC_cohort_multires_harmony/data/scVKMYC_cohort_multires_harmony_theta05_seurat.rds"

fileName <- "scVKMYC_cohort_multires_harmony_theta05_seurat"

seurat <- readRDS(seuratObj)


###########################################################
### UMAP/TSNE PLOT - colour by cluster/sample
###########################################################

meta <- seurat@meta.data
meta$Sample_ID <- factor(as.character(meta$Sample_ID),
                         levels = c("Cont1", "Cont2", "Cont3",
                                     "ND1", "ND2", "ND3", "ND4", "ND5",
                                     "AMG1", "AMG2", "AMG3",
                                     "MM1", "MM2", "MM3","MM4","MM5","MM6","MM7")) 

## -------- By sample ------- ##


ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Sample_ID")) +
  geom_point(size = 1, alpha = 0.8) + theme_classic() +
  theme_bw() + theme(axis.text.x = element_blank(), 
                     axis.text.y = element_blank(), 
                     axis.ticks = element_blank(),
                     panel.border = element_rect(linetype = "solid",
                                                fill = NA, size = 1),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank()) +
  labs(x="tSNE_1", y="tSNE_2") +
  ggtitle("B Lineage Cells (n=17,504) \nColoured by Sample") 

```


### Figure 1C. SDC1/CD19 Expression

```{r}

###########################################################
### Load input data
###########################################################


seuratObj <- "/Volumes/Samsung_T5/OneDrive - UHN/Documents/scVKMYC/analysis/cohort/lineages/scVKMYC_cohort_lineages/B_PC/scVKMYC_cohort_multires_harmony/data/scVKMYC_cohort_multires_harmony_theta05_seurat.rds"

fileName <- "scVKMYC_cohort_multires_harmony_theta05_seurat"

seurat <- readRDS(seuratObj)

###########################################################
### CD19/CD138 Co-expression Plots
###########################################################


#co-expression plots 

FeaturePlot(object = seurat, c("Cd19", "Sdc1"), reduction = "tsne",
    cols = c("grey", "blue", "red"), 
    blend = TRUE) +
  theme_bw() + theme(axis.text.x = element_blank(), 
                         axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid",
                                                    fill = NA, size = 1),
                        # panel.border = element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.position = "none",
                        title=element_text(face="italic"))

```



### Figure 1C. Plasma cell score

```{r}

###########################################################
### Load input data
###########################################################


seuratObj <- "/Volumes/Samsung_T5/OneDrive - UHN/Documents/scVKMYC/analysis/cohort/lineages/scVKMYC_cohort_lineages/B_PC/scVKMYC_cohort_multires_harmony/data/scVKMYC_cohort_multires_harmony_theta05_seurat.rds"

fileName <- "scVKMYC_cohort_multires_harmony_theta05_seurat"

seurat <- readRDS(seuratObj)


###########################################################
### Plasma cell gene set score 
###########################################################

dat <- seurat

runSignature <- "HAY_BONE_MARROW_PLASMA_CELL_MOUSE"

geneSignatures <- "~/OneDrive - UHN/Documents/10x_Experiments/MASTER_genesets.csv"
sigs <- read.csv(geneSignatures, stringsAsFactors = F, blank.lines.skip = T, na.strings = "")
runSig <- sigs[2:nrow(sigs), ,drop = F]
runSig <- runSig[, runSignature]
runSig <- runSig[!is.na(runSig)]

#add title and length of signature to name 
sigs <- list(runSig)
names(sigs) <- runSignature

dat <- AddModuleScore(dat,
                     features = sigs,
                     ctrl = 25,
                     name = make.names(names(sigs))
                    )
start <- ncol(dat@meta.data)-length(sigs) + 1
end <- ncol(dat@meta.data)
colnames(dat@meta.data)[start:end] <- paste0(names(sigs), "_ModuleScore")


FeaturePlot(object = dat, c("HAY_BONE_MARROW_PLASMA_CELL_MOUSE_ModuleScore"), reduction = "tsne") +
  theme_bw() + theme(axis.text.x = element_blank(), 
                         axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid",
                                                    fill = NA, size = 1),
                        # panel.border = element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.position = "none",
                        title=element_text(face="italic"))



```


### Figure 1G. MYC transgene expression

```{r}

seurat.file <- paste0("./data/", "scVKMYC_TSNErefined_PConly", "_seurat.rds")
dat <- readRDS(seurat.file)

FeaturePlot(object = dat, c("MYC"), reduction = "tsne") +
  theme_bw() + theme(axis.text.x = element_blank(), 
                         axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                     axis.title = element_blank(),
                        panel.border = element_blank(),
                        # panel.border = element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.position = "none"
                        )




```


