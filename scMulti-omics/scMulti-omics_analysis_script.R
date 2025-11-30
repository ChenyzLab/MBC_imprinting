# Preprocessing and analysis of scMulti-omics data, some of the code and analytical strategies used in this script were adapted from https://doi.org/10.1038/s41590-023-01497-y

# This script is divided into six parts:
# Part 1: Loading R packages and functions, and defining the colors used in this paper.
# Part 2 and 3: Preprocessing of two batch datasets, including demultiplexing, quality control (QC), non-B cell filtering, and defining antigen+ B cells.
# Part 4: Dataset integration, WNN analysis, and cell type annotation.
# Part 5: Identifying differentially expressed genes and performing pathway analysis.
# Part 6: BCR repertoire analysis

#### 01. Load packages and functions, define colors ####
library(seqinr)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(bruceR)
library(reshape2)
library(data.table)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(graphics)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(colorspace)
library(ggridges)
library(scales)
library(harmony)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(readxl)
library(Biobase)
library(gridExtra)
library(ggExtra)
library(rstatix)
library(GSEABase)


load("My_func.RData")
BC_col <- c("rNAV" = "#4F99C9","aNAV" = "#A6D0E6","DN2" = "#FFC0CB","DN1" = "#c0937e","uMBC" = "#A8D3A0","sMBC" = "#8582BD","aMBC" = "#FA8072","PB" = "#8B4789")
Antigen <- c("S2" = "#A6D0E6","BF7 RBD" = "#b71f49","Cross RBD" = "#f58d60","Ancestral RBD" =  "#EEDC82",
             "BF7 NTD" = "#FA8090", "Cross NTD" = "#f9ab60","Ancestral NTD"= "#CDBE70", "Ag-" = "grey95")
HI = "#839FBF"
HVI = "#074080"
SLEI = "#DFBFDF"
SLEVI = "#800080"
BF7 = "#88C0C0"
Cross = "#E4E4E4"
Ancestral = "grey95"
Early = "#D4B81C"
Late = "#396060"

#### 02. Batch1 preprocessing ####
#### 02.1 Batch1: loading data, filtering, Hashing demultiplexing, add metadata information ####
cell <- Read10X("./cellranger/batch1/sample_filtered_feature_bc_matrix")
raw <- CreateSeuratObject(counts = cell$`Gene Expression`)
rownames(cell$`Antibody Capture`)
#[1] "CD21"     "FCRL5.1"  "CD71"     "CXCR5.1"  "IgG"      "IgM"      "CD11C"    "CD27.1"   "IgD"      "CD38.1"  
#[11] "Hashtag1" "Hashtag2" "Hashtag3" "Hashtag4" "Hashtag5" "Hashtag6" "Hashtag7" "Hashtag8" "APC1"     "APC2.1"  
#[21] "APC3"     "APC4"     "APC5"     "APC6"     "APC7" 

# Creat assays for Baiting, Protein and hashing
Baiting_assay <- CreateAssayObject(counts = cell$`Antibody Capture`[c(19:25),])
Protein_assay <- CreateAssayObject(counts = cell$`Antibody Capture`[c(1:10),])
Hashing_assay <- CreateAssayObject(counts = cell$`Antibody Capture`[c(11:18),])

# Add assays to the previously created Seurat object
raw[["Baiting"]] <- Baiting_assay
raw[["Protein"]] <- Protein_assay
raw[["Hashing"]] <- Hashing_assay

remove(Baiting_assay,Hashing_assay,Protein_assay,cell)

# Demultiplexing the HTO data
raw <- NormalizeData(raw, assay = "Hashing", normalization.method = "CLR", margin = 2)
raw <- HTODemux(raw, assay = "Hashing", positive.quantile = 0.99)
RidgePlot(raw, assay = "Hashing", features = rownames(raw[["Hashing"]]), ncol = 4)+center.title()
table(raw$Hashing_classification.global)

# Filter out cells that are not singlets
singlet <- subset(raw, Hashing_classification.global == "Singlet")

# QC
Ig_gene <- grep("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]",rownames(singlet),value = T)

singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
singlet[["percent.Ig"]] <- PercentageFeatureSet(singlet, features = Ig_gene)
singlet[["percent.ribo"]] <- PercentageFeatureSet(singlet, pattern = "^RP[SL]")

Idents(singlet) <- "orig.ident"
VlnPlot(singlet, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3,pt.size = 0)

quantile(singlet@meta.data$nFeature_RNA,probs = seq(0, 1, 1/100))
quantile(singlet@meta.data$percent.mt,probs = seq(0, 1, 1/100))
singlet <- subset(singlet, subset = nFeature_RNA > 200 & percent.mt < 10)

singlet <- RenameCells(singlet, new.names = paste0(colnames(singlet), "_batch1"))
singlet$cell_id <- rownames(singlet@meta.data)
rm(raw)

# Remove Ig genes
singlet@assays$RNA <- subset(singlet@assays$RNA, features=setdiff(rownames(singlet@assays$RNA), Ig_gene)) 

# Add timepoint and sample information
singlet$sampleid <- ""
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag1"] <- "H002-1"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag2"] <- "H002-2"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag3"] <- "H050-1"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag4"] <- "H050-2"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag5"] <- "SLE_ZJJ"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag6"] <- "SLE_WL"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag7"] <- "CR-11"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag8"] <- "SLE_LYL"

singlet$sampleinfo <- ""
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("H002-1","H002-2")] <- "HVI"
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("H050-1","H050-2")] <- "HI"
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("SLE_ZJJ","SLE_WL")] <- "SLEVI"
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("CR-11","SLE_LYL")] <- "SLEI"

singlet$timepoint <- ""
singlet$timepoint[singlet@meta.data$sampleid %in% c("H002-1","H050-1","SLE_ZJJ","CR-11")] <- "early"
singlet$timepoint[singlet@meta.data$sampleid %in% c("H002-2","H050-2","SLE_WL","SLE_LYL")] <- "late"

table(singlet$sampleid)
table(singlet$sampleinfo)
table(singlet$timepoint)

# We remove SLE_LYL because the infection timepoint is uncetain
singlet <- subset(singlet,sampleid != "SLE_LYL")

#### 02.2 Batch1: Reduce dimension and Cluster ####
# processing step: Normalization, scale, reduce dimension and cluster
DefaultAssay(singlet) <- "RNA"
singlet <- NormalizeData(singlet)
singlet <- FindVariableFeatures(singlet)
singlet <- ScaleData(singlet)
singlet <- RunPCA(singlet)
ElbowPlot(singlet,50)
singlet <- FindNeighbors(singlet, dims = 1:30)
singlet <- RunUMAP(singlet, dims = 1:30)
singlet <- FindClusters(singlet, resolution = 0.5)

marker <- c("ITGAM","CD3D", "CD3E", "CD3G","CD4","CD8A","CD14","NKG7","NCAM1","HBB","CD1C","FCER1A","IL7R","RAG1","RAG2","MME","CD79A","CD79B","MS4A1","CD19","CR2","CD27","CD38","TFRC","CXCR5","ITGAX","FCRL5","SDC1","PRDM1")
My_DotPlot(singlet,features = c(marker),group_by = "seurat_clusters",dot_scale = 6,coord_flip = FALSE)

# remove c11,14 for doublets and rerun the processing step above
singlet <- subset(singlet, seurat_clusters %in% c(11,14),invert = T)
singlet <- NormalizeData(singlet)
singlet <- FindVariableFeatures(singlet)
singlet <- ScaleData(singlet)
singlet <- RunPCA(singlet)
ElbowPlot(singlet,50)
singlet <- FindNeighbors(singlet, dims = 1:30)
singlet <- RunUMAP(singlet, dims = 1:30)
singlet <- FindClusters(singlet, resolution = 0.5)

#Then, annotate main cell type with new cluster
singlet$mainType <- "BC"
singlet$mainType[singlet$seurat_clusters %in% c(7,12)] <- "Mono_DC"
singlet$mainType[singlet$seurat_clusters %in% c(9,11)] <- "T_NK"
singlet$cell_id <- rownames(singlet@meta.data)
table(singlet$mainType)

# filter nonB cells and subset B cells for downstream analysis
BC <- subset(singlet, mainType == "BC") #10181 cells
BC <- NormalizeData(BC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(singlet,50)
BC <- FindNeighbors(BC, dims = 1:30)
BC <- RunUMAP(BC, dims = 1:30)
BC <- FindClusters(BC, resolution = 0.5)
My_DotPlot(BC,features = c(marker),group_by = "seurat_clusters",dot_scale = 6,coord_flip = FALSE)

#### 02.3 Batch1: Processing of baiting counts ####
# Pre-process baiting counts
Baiting_df <- as.data.frame(t(as.data.frame(BC@assays$Baiting@counts)))
colnames(Baiting_df) <- c("Ancestral_RBD", "BF7_RBD","Ancestral_NTD","BF7_NTD","S2","HA","Blank")

# Density plots
melted_Baiting_df <- reshape2::melt(Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 6)+
  scale_x_continuous(trans='log2',breaks = c(0,1,2,3,4,5,10,20,40,100,200,500,1000,3000))+
  labs(title = "Density plots before negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())

# Correct by using the negative control
Baiting_df$cor_BF7_RBD <- Baiting_df$BF7_RBD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_BF7_NTD <- Baiting_df$BF7_NTD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_Ancestral_RBD <- Baiting_df$Ancestral_RBD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_Ancestral_NTD <- Baiting_df$Ancestral_NTD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_S2 <- Baiting_df$S2 - pmax(Baiting_df$Blank,Baiting_df$HA)

# Density plots after negative control subtraction
melted_Baiting_df <- Baiting_df[c(12,8:11)]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
p <- ggplot(melted_Baiting_df[melted_Baiting_df$value>0,], aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+ 
  scale_fill_manual(values = c("#A6D0E6","#b71f49","#FA8090","#EEDC82","#CDBE70")) +
  scale_x_continuous(trans='log2',breaks = c(1,2,5,8,14,19,40,100,400))+
  labs(title = "Density plots define Antigen-reactive cutoffs", x="Reads", y = "") +
  theme_ridges(font_size = 25, grid = TRUE) + 
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+NoLegend()
ggsave("./for_publication/figure/LIBRAseq/FigS4C.Batch1_cutoff.pdf",p,width = 8,height = 6)

# visualize density plot
melted_Baiting_df <- Baiting_df[8:12]
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_S2 > 0, "cor_S2"])#cutoff:5.551247
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_BF7_RBD > 0, "cor_BF7_RBD"]) #cutoff:5.585893
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_BF7_NTD > 0, "cor_BF7_NTD"]) #cutoff:19.807030
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_Ancestral_RBD > 0, "cor_Ancestral_RBD"])#cutoff:14.153433
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_Ancestral_NTD > 0, "cor_Ancestral_NTD"])#cutoff:8.646832

# Classify binders
Baiting_df$cor_BF7_RBD[Baiting_df$cor_BF7_RBD<1] <- 0
Baiting_df$cor_BF7_NTD[Baiting_df$cor_BF7_NTD<1] <- 0
Baiting_df$cor_Ancestral_RBD[Baiting_df$cor_Ancestral_RBD<1] <- 0
Baiting_df$cor_Ancestral_NTD[Baiting_df$cor_Ancestral_NTD<1] <- 0
Baiting_df$cor_S2[Baiting_df$cor_S2<1] <- 0

Baiting_df$cor_BF7_RBD.classification <- "Negative"
Baiting_df$cor_BF7_NTD.classification <- "Negative"
Baiting_df$cor_Ancestral_RBD.classification <- "Negative"
Baiting_df$cor_Ancestral_NTD.classification <- "Negative"
Baiting_df$cor_S2.classification <- "Negative"

# Determine the cut-off values for antigen specificity based on the bimodal distribution
Baiting_df$cor_BF7_RBD.classification[Baiting_df$cor_BF7_RBD>=5] <- "Positive"
Baiting_df$cor_BF7_NTD.classification[Baiting_df$cor_BF7_NTD>=19] <- "Positive"
Baiting_df$cor_Ancestral_RBD.classification[Baiting_df$cor_Ancestral_RBD>=14] <- "Positive"
Baiting_df$cor_Ancestral_NTD.classification[Baiting_df$cor_Ancestral_NTD>8] <- "Positive"
Baiting_df$cor_S2.classification[Baiting_df$cor_S2>=5] <- "Positive"

# Adding a column about whether a cell is positive for antigens
Baiting_df$bait.positive <- "no"

Baiting_df$bait.positive[Baiting_df$cor_S2.classification == "Positive" |
                           Baiting_df$cor_BF7_NTD.classification == "Positive" |
                           Baiting_df$cor_Ancestral_RBD.classification == "Positive" |
                           Baiting_df$cor_BF7_RBD.classification == "Positive" |
                           Baiting_df$cor_Ancestral_NTD.classification == "Positive"] <- "yes"
table(Baiting_df$bait.positive)

# Setting corrected baiting counts to NA if classified as binding negative
Baiting_df_forcal <- Baiting_df
Baiting_df_forcal$cor_BF7_RBD[Baiting_df_forcal$cor_BF7_RBD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_BF7_NTD[Baiting_df_forcal$cor_BF7_NTD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_Ancestral_RBD[Baiting_df_forcal$cor_Ancestral_RBD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_Ancestral_NTD[Baiting_df_forcal$cor_Ancestral_NTD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_S2[Baiting_df_forcal$cor_S2.classification == "Negative"] <- NA

# Density plots after cutoff correction
melted_Baiting_df_forcal <- Baiting_df_forcal[8:12]
melted_Baiting_df_forcal <- reshape2::melt(melted_Baiting_df_forcal)
melted_Baiting_df_forcal <- na.omit(melted_Baiting_df_forcal)

ggplot(melted_Baiting_df_forcal, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 1)+
  scale_fill_manual(values = c("#b71f49","#FA8090","#EEDC82","#CDBE70","#A6D0E6")) +
  scale_x_continuous(trans='log2',breaks = c(1,2,3,4,5,10,20,40,100,200,500))+
  labs(title = "Density plots after cutoff correction") +
  labs(title = "Batch1: Density plots after cutoff correction", x="Reads", y = "") +
  theme_ridges(font_size = 25, grid = TRUE) + 
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+NoLegend()
remove(melted_Baiting_df_forcal)

# Normalize and scale of corrected baiting assay
Baiting_df_forcal <- as.data.frame(t(Baiting_df_forcal[,8:12]))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df_forcal)
BC[["Cor_Baiting"]] <- Cor_Baiting_assay

remove(Cor_Baiting_assay)

BC <- NormalizeData(BC, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1) #normalized across antigen
normalized_across_features <- as.data.frame(t(BC@assays$Cor_Baiting@data))
colnames(normalized_across_features) <- c("cor_BF7_RBD","cor_BF7_NTD","cor_Ancestral_RBD","cor_Ancestral_NTD","cor_S2")
remove(Baiting_df_forcal)
BC@meta.data$nCount_Cor_Baiting <- NULL
BC@meta.data$nFeature_Cor_Baiting <- NULL
normalized_across_features <- myscaling(normalized_across_features)  

normalized_across_features.melted <- reshape2::melt(normalized_across_features)
normalized_across_features.melted <- na.omit(normalized_across_features.melted)
normalized_across_features.melted <- normalized_across_features.melted[normalized_across_features.melted$value>0,]
ggplot(normalized_across_features.melted, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  labs(title = "Density plots after scaling") + 
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(normalized_across_features.melted)

colnames(normalized_across_features) <- c("BF7_RBD_score","BF7_NTD_score","Ancestral_RBD_score","Ancestral_NTD_score","S2_score")

Baiting_df_final <- merge(Baiting_df,normalized_across_features,by = 0)
rownames(Baiting_df_final) <- Baiting_df_final$Row.names
Baiting_df_final$Row.names <- NULL

# These are the final LIBRA scores, adding these scores to the Seurat object
BC@meta.data <- merge(BC@meta.data,Baiting_df_final, by=0)
rownames(BC@meta.data) <- BC$Row.names
BC$Row.names <- NULL
remove(normalized_across_features,Baiting_df_final,Baiting_df)

# Classify antigens
# If one cell is reactive to both S2, Ancestral(WH) and BF.7, only if the score of S2 is higher than Ancestral(WH) and BF。7, it is classified as S2 
BC$antigen <- "Ag-"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Positive" & BC$cor_BF7_RBD.classification == "Negative" & BC$cor_S2.classification == "Negative"] <- "Ancestral RBD"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Negative" & BC$cor_BF7_RBD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "BF7 RBD"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Positive" & BC$cor_BF7_RBD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "Cross RBD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Positive" & BC$cor_BF7_NTD.classification == "Negative" & BC$cor_S2.classification == "Negative"] <- "Ancestral NTD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Negative" & BC$cor_BF7_NTD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "BF7 NTD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Positive" & BC$cor_BF7_NTD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "Cross NTD"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Negative" & BC$cor_BF7_RBD.classification == "Negative" & 
             BC$cor_Ancestral_NTD.classification == "Negative" & BC$cor_BF7_NTD.classification == "Negative" & BC$cor_S2.classification == "Positive"] <- "S2"
BC$antigen[BC$cor_BF7_RBD.classification == "Positive" & BC$cor_Ancestral_RBD.classification == "Positive" & 
             BC$cor_S2.classification == "Positive" & BC$BF7_RBD_score > BC$S2_score] <- "Cross RBD"

table(BC$bait.positive,BC$antigen)

rm(melted_Baiting_df)

#### 02.4 Batch1: add BCR info, classify isotype ####
BCR_batch1 <- read.table("./cellranger/batch1/vdj/filtered_contig.tsv",header = T,sep = "\t")
BCR_batch1$cell_id <- paste0(sapply(strsplit(BCR_batch1$sequence_id,"_"),'[',1),"_batch1")
BCR_anno <- read.csv("./cellranger/batch1/vdj/filtered_contig_annotations.csv")
BCR_anno <- BCR_anno[,c(3,27,28)]
colnames(BCR_anno) <- c("sequence_id","reads","umis")

BCR <- merge(BCR_batch1,BCR_anno,by = "sequence_id",all.x = T)
rm(BCR_batch1,BCR_anno)

BCR <- BCR %>%
  mutate(chain = ifelse(locus == 'IGH', 'Heavy', 'Light')) %>%
  group_by(cell_id) %>%
  mutate(
    type = case_when(
      sum(chain == "Heavy") >= 2 ~ "double",
      sum(chain == "Light") >= 2 ~ "double",
      sum(chain == "Heavy") == 0 & sum(chain == "Light") < 2 ~ "noHeavy",
      sum(chain == "Light") == 0 & sum(chain == "Heavy") < 2 ~ "noLight",
      TRUE ~ "single"
    )
  ) %>%
  ungroup()
table(BCR$type)

# double noHeavy noLight  single 
# 2951    1694     389    20434

# some B cells have more than one heavy or light chain, I chose the dominant one.
BCR <- BCR %>%
  mutate(chain = ifelse(locus == 'IGH', 'Heavy', 'Light')) %>%
  group_by(cell_id, chain) %>%
  arrange(cell_id, desc(productive), desc(umis), desc(reads)) %>%
  slice_head(n=1)

BCR_H <- BCR[BCR$chain == "Heavy",]
BCR_H$SHM_count <- round((BCR_H$v_sequence_end - BCR_H$v_sequence_start +1) * (1 - (BCR_H$v_identity)/100))
BCR_H$SHM_freq <- 1 - BCR_H$v_identity/100
BCR_H$cdr3_length <- BCR_H$cdr3_end - BCR_H$cdr3_start + 1
colnames(BCR_H)[colnames(BCR_H) != "cell_id"] <- paste0("BCR_Heavy_", colnames(BCR_H)[colnames(BCR_H) != "cell_id"])

BCR_L <- BCR[BCR$chain == "Light",]
BCR_L$SHM_count <- round((BCR_L$v_sequence_end - BCR_L$v_sequence_start +1) * (1 - (BCR_L$v_identity)/100))
BCR_L$SHM_freq <- 1 - BCR_L$v_identity/100
BCR_L$cdr3_length <- BCR_L$cdr3_end - BCR_L$cdr3_start + 1
colnames(BCR_L)[colnames(BCR_L) != "cell_id"] <- paste0("BCR_Light_", colnames(BCR_L)[colnames(BCR_L) != "cell_id"])
BCR_merge_batch1 <- merge(BCR_H ,BCR_L, by='cell_id', all.x = T)

metadata <- BC@meta.data
metadata <- merge(metadata,BCR_merge_batch1,by = "cell_id",all.x = T)
rownames(metadata) <- metadata$cell_id
order <- colnames(BC)
metadata <- metadata[order,]

BC@meta.data <- metadata

BC <- add_BCR_info(BC, classify_col = "BCR_Heavy_c_call", new_col = "sub_Isotype")
table(BC$sub_Isotype)

BC@meta.data$Isotype <- "unknown"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype %in% c("IGHG1", "IGHG2","IGHG3","IGHG4","IGHG4A")] <- "IGHG"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype %in% c("IGHA1","IGHA2")] <- "IGHA"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype == "IGHD"] <- "IGHD"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype == "IGHM"] <- "IGHM"
table(BC$Isotype)

BC$RNA_snn_res.0.5 <- NULL
saveRDS(BC,"./for_publication/result/BC_batch1.rds")
rm(BCR,BCR_H,BCR_L,BCR_merge_batch1,metadata,singlet,BC)

#### 03. Batch2 preprocesing ####
#### 03.1 Batch2: loading data, filtering, Hashing demultiplexing, add metadata information ####
cell <- Read10X("./cellranger/batch2/sample_filtered_feature_bc_matrix")
raw <- CreateSeuratObject(counts = cell$`Gene Expression`)
rownames(cell$`Antibody Capture`)
#[1] "CD21"     "FCRL5.1"  "CD71"     "CXCR5.1"  "IgG"      "IgM"      "CD11C"    "CD27.1"   "IgD"      "CD38.1"   "Hashtag1"
#[12] "Hashtag2" "Hashtag3" "Hashtag4" "Hashtag5" "Hashtag6" "Hashtag7" "Hashtag8" "APC1"     "APC2.1"   "APC3"     "APC4"    
#[23] "APC5"     "APC6"     "APC7"

# Creat assays for Baiting, Protein and hashing
Baiting_assay <- CreateAssayObject(counts = cell$`Antibody Capture`[c(19:25),])
Protein_assay <- CreateAssayObject(counts = cell$`Antibody Capture`[c(1:10),])
Hashing_assay <- CreateAssayObject(counts = cell$`Antibody Capture`[c(11:18),])

# Add assays to the previously created Seurat object
raw[["Baiting"]] <- Baiting_assay
raw[["Protein"]] <- Protein_assay
raw[["Hashing"]] <- Hashing_assay

remove(Baiting_assay,Hashing_assay,Protein_assay,cell)

# Demultiplex the HTO data
raw <- NormalizeData(raw, assay = "Hashing", normalization.method = "CLR", margin = 2) #normalize across cells
raw <- HTODemux(raw, assay = "Hashing", positive.quantile = 0.99)
RidgePlot(raw, assay = "Hashing", features = rownames(raw[["Hashing"]]), ncol = 4)+center.title()
table(raw$Hashing_classification.global)
#Doublet Negative  Singlet 
#5900     2876    29449
table(raw$Hashing_classification.global,raw$Hashing_maxID)

# Filter out cells that are not singlets
singlet <- subset(raw, Hashing_classification.global=="Singlet")

# QC
Ig_gene <- grep("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]",rownames(singlet),value = T)

singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
singlet[["percent.Ig"]] <- PercentageFeatureSet(singlet, features = Ig_gene)
singlet[["percent.ribo"]] <- PercentageFeatureSet(singlet, pattern = "^RP[SL]")

Idents(singlet) <- "orig.ident"
VlnPlot(singlet, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3,pt.size = 0)

quantile(singlet@meta.data$nFeature_RNA,probs = seq(0, 1, 1/100))
quantile(singlet@meta.data$percent.mt,probs = seq(0, 1, 1/100))
singlet <- subset(singlet, subset = nFeature_RNA > 200 & percent.mt < 10) #29313 cells

singlet <- RenameCells(singlet, new.names = paste0(colnames(singlet), "_batch2"))
singlet$cell_id <- rownames(singlet@meta.data)

rm(raw)

# Remove Ig genes
singlet@assays$RNA <- subset(singlet@assays$RNA, features=setdiff(rownames(singlet@assays$RNA), Ig_gene)) 

# Add timepoint and sample information
singlet$sampleid <- ""
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag1"] <- "H033-1"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag2"] <- "H033-2"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag3"] <- "H052-1"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag4"] <- "H052-2"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag5"] <- "R014-1"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag6"] <- "CR-3"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag7"] <- "R017-1"
singlet$sampleid[singlet@meta.data$hash.ID == "Hashtag8"] <- "R017-2"

singlet$sampleinfo <- ""
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("H033-1","H033-2")] <- "HVI"
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("H052-1","H052-2")] <- "HI"
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("R014-1","CR-3")] <- "SLEVI"
singlet$sampleinfo[singlet@meta.data$sampleid %in% c("R017-1","R017-2")] <- "SLEI"

singlet$timepoint <- ""
singlet$timepoint[singlet@meta.data$sampleid %in% c("H033-1","H052-1","R014-1","R017-1")] <- "early"
singlet$timepoint[singlet@meta.data$sampleid %in% c("H033-2","H052-2","CR-3","R017-2")] <- "late"

table(singlet$sampleid)

#### 03.2 Batch2: Reduce dimension and Cluster ####
# processing step: Normalization, scale, reduce dimension and cluster
DefaultAssay(singlet) <- "RNA"
singlet <- NormalizeData(singlet)
singlet <- FindVariableFeatures(singlet)
singlet <- ScaleData(singlet)
singlet <- RunPCA(singlet)
ElbowPlot(singlet,50)
singlet <- FindNeighbors(singlet, dims = 1:30)
singlet <- RunUMAP(singlet, dims = 1:30)
singlet <- FindClusters(singlet, resolution = 0.5)

DimPlot(singlet, reduction = 'umap', label = T, label.size = 5, pt.size = 0.8,group.by = c("seurat_clusters","sampleid")) 

marker <- c("ITGAM","CD3D", "CD3E", "CD3G","CD4","CD8A","CD14","NKG7","NCAM1","HBB","CD1C","FCER1A","IL7R","RAG1","RAG2","MME","CD79A","CD79B","MS4A1","CD19","CR2","CD27","CD38","TFRC","CXCR5","ITGAX","FCRL5","SDC1","PRDM1") #Itgam is Cd11b(Myeloid Cells)
My_DotPlot(singlet,features = c(marker),group_by = "seurat_clusters",dot_scale = 6,coord_flip = FALSE)

# remove c10,18 for doublets and rerun the processing step above
singlet <- subset(singlet, seurat_clusters %in% c(10,18),invert = T)
singlet <- NormalizeData(singlet)
singlet <- FindVariableFeatures(singlet)
singlet <- ScaleData(singlet)
singlet <- RunPCA(singlet)
ElbowPlot(singlet,50)
singlet <- FindNeighbors(singlet, dims = 1:30)
singlet <- RunUMAP(singlet, dims = 1:30)
singlet <- FindClusters(singlet, resolution = 0.5)
#Then, annotate main cell type with new cluster
singlet$mainType <- "BC"
singlet$mainType[singlet$seurat_clusters %in% c(9,11)] <- "Mono_DC"
singlet$mainType[singlet$seurat_clusters %in% c(6,7,10,12,15)] <- "T_NK"
singlet$cell_id <- rownames(singlet@meta.data)
table(singlet$mainType)

# subset BC
BC <- subset(singlet, mainType == "BC")
BC <- NormalizeData(BC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(BC,50)
BC <- FindNeighbors(BC, dims = 1:30)
BC <- RunUMAP(BC, dims = 1:30)
BC <- FindClusters(BC, resolution = 0.5)

DimPlot(BC, reduction = 'umap', label = T, label.size = 4, pt.size = 0.8,group.by = c("seurat_clusters","sampleid"))

My_DotPlot(BC, features = c(marker), group_by = "seurat_clusters",dot_scale = 6,coord_flip = FALSE)
VlnPlot(BC, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.Ig"), ncol = 2,pt.size = 0,group.by = "seurat_clusters")

# remove c9 for doublets and rerun the processing pipline above
BC <- subset(BC, seurat_clusters != 9) #21309 cells

#### 03.3 Batch2: Processing of baiting counts ####
# Pre-processing baiting counts
Baiting_df <- as.data.frame(t(as.data.frame(BC@assays$Baiting@counts)))
colnames(Baiting_df) <- c("Ancestral_RBD", "BF7_RBD","Ancestral_NTD","BF7_NTD","S2","HA","Blank")

# Density plots
melted_Baiting_df <- reshape2::melt(Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 6)+
  scale_x_continuous(trans='log2',breaks = c(0,1,2,3,4,5,10,20,40,100,200,500,1000,3000))+
  labs(title = "Density plots before negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())

# Correcting by using the negative control
Baiting_df$cor_BF7_RBD <- Baiting_df$BF7_RBD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_BF7_NTD <- Baiting_df$BF7_NTD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_Ancestral_RBD <- Baiting_df$Ancestral_RBD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_Ancestral_NTD <- Baiting_df$Ancestral_NTD - pmax(Baiting_df$Blank,Baiting_df$HA)
Baiting_df$cor_S2 <- Baiting_df$S2 - pmax(Baiting_df$Blank,Baiting_df$HA)

# Density plots after negative control subtraction
melted_Baiting_df <- Baiting_df[c(12,8:11)]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
p <- ggplot(melted_Baiting_df[melted_Baiting_df$value>1,], aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+ 
  scale_fill_manual(values = c("#A6D0E6","#b71f49","#FA8090","#EEDC82","#CDBE70")) +
  scale_x_continuous(trans='log2',breaks = c(1,2,5,6,16,40,100,400))+
  labs(title = "Batch2: Density plots after background correction", x="Reads", y = "") +
  theme_ridges(font_size = 25, grid = TRUE) + 
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+NoLegend()

# visualize density plot
melted_Baiting_df <- Baiting_df[8:12]

cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_BF7_RBD > 0, "cor_BF7_RBD"]) #cutoff:5.751285
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_BF7_NTD > 0, "cor_BF7_NTD"]) #cutoff:6.783424
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_Ancestral_RBD > 0, "cor_Ancestral_RBD"]) #cutoff:40.142773
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_Ancestral_NTD > 0, "cor_Ancestral_NTD"]) #cutoff:16.114073
cal_density_cutoff(melted_Baiting_df[melted_Baiting_df$cor_S2 > 0, "cor_S2"])#cutoff:6.969090

# Classify binders
Baiting_df$cor_BF7_RBD[Baiting_df$cor_BF7_RBD<1] <- 0
Baiting_df$cor_BF7_NTD[Baiting_df$cor_BF7_NTD<1] <- 0
Baiting_df$cor_Ancestral_RBD[Baiting_df$cor_Ancestral_RBD<1] <- 0
Baiting_df$cor_Ancestral_NTD[Baiting_df$cor_Ancestral_NTD<1] <- 0
Baiting_df$cor_S2[Baiting_df$cor_S2<1] <- 0

Baiting_df$cor_BF7_RBD.classification <- "Negative"
Baiting_df$cor_BF7_NTD.classification <- "Negative"
Baiting_df$cor_Ancestral_RBD.classification <- "Negative"
Baiting_df$cor_Ancestral_NTD.classification <- "Negative"
Baiting_df$cor_S2.classification <- "Negative"

# Determine the cut-off values for antigen specificity based on the bimodal distribution
Baiting_df$cor_BF7_RBD.classification[Baiting_df$cor_BF7_RBD>=5] <- "Positive"
Baiting_df$cor_BF7_NTD.classification[Baiting_df$cor_BF7_NTD>=6] <- "Positive"
Baiting_df$cor_Ancestral_RBD.classification[Baiting_df$cor_Ancestral_RBD>=40] <- "Positive"
Baiting_df$cor_Ancestral_NTD.classification[Baiting_df$cor_Ancestral_NTD>=16] <- "Positive"
Baiting_df$cor_S2.classification[Baiting_df$cor_S2>=6] <- "Positive"

# Adding a column about whether a cell is positive for antigens
Baiting_df$bait.positive <- "no"

Baiting_df$bait.positive[Baiting_df$cor_S2.classification == "Positive" |
                           Baiting_df$cor_BF7_NTD.classification == "Positive" |
                           Baiting_df$cor_Ancestral_RBD.classification == "Positive" |
                           Baiting_df$cor_BF7_RBD.classification == "Positive" |
                           Baiting_df$cor_Ancestral_NTD.classification == "Positive"] <- "yes"
table(Baiting_df$bait.positive)

# Setting corrected baiting counts to NA if classified as binding negative
Baiting_df_forcal <- Baiting_df
Baiting_df_forcal$cor_BF7_RBD[Baiting_df_forcal$cor_BF7_RBD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_BF7_NTD[Baiting_df_forcal$cor_BF7_NTD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_Ancestral_RBD[Baiting_df_forcal$cor_Ancestral_RBD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_Ancestral_NTD[Baiting_df_forcal$cor_Ancestral_NTD.classification == "Negative"] <- NA
Baiting_df_forcal$cor_S2[Baiting_df_forcal$cor_S2.classification == "Negative"] <- NA


# Density plots after cutoff correction
melted_Baiting_df_forcal <- Baiting_df_forcal[8:12]
melted_Baiting_df_forcal <- reshape2::melt(melted_Baiting_df_forcal)
melted_Baiting_df_forcal <- na.omit(melted_Baiting_df_forcal)

ggplot(melted_Baiting_df_forcal, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 1)+
  scale_fill_manual(values = c("#b71f49","#FA8090","#EEDC82","#CDBE70","#A6D0E6")) +
  scale_x_continuous(trans='log2',breaks = c(1,2,3,4,5,10,20,40,100,200,500))+
  labs(title = "Density plots after cutoff correction") +
  labs(title = "Batch2: Density plots after cutoff correction", x="Reads", y = "") +
  theme_ridges(font_size = 25, grid = TRUE) + 
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+NoLegend()
remove(melted_Baiting_df_forcal,dens)

# Normalize and scale of corrected baiting assay
Baiting_df_forcal <- as.data.frame(t(Baiting_df_forcal[,8:12]))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df_forcal)
BC[["Cor_Baiting"]] <- Cor_Baiting_assay

remove(Cor_Baiting_assay)

BC <- NormalizeData(BC, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1) #normalized across antigen
normalized_across_features <- as.data.frame(t(BC@assays$Cor_Baiting@data))
colnames(normalized_across_features) <- c("cor_BF7_RBD","cor_BF7_NTD","cor_Ancestral_RBD","cor_Ancestral_NTD","cor_S2")
remove(Baiting_df_forcal)
BC@meta.data$nCount_Cor_Baiting <- NULL
BC@meta.data$nFeature_Cor_Baiting <- NULL
normalized_across_features <- myscaling(normalized_across_features)  

normalized_across_features.melted <- reshape2::melt(normalized_across_features)
normalized_across_features.melted <- na.omit(normalized_across_features.melted)
normalized_across_features.melted <- normalized_across_features.melted[normalized_across_features.melted$value>0,]
ggplot(normalized_across_features.melted, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  labs(title = "Density plots after scaling") + 
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(normalized_across_features.melted)

colnames(normalized_across_features) <- c("BF7_RBD_score","BF7_NTD_score","Ancestral_RBD_score","Ancestral_NTD_score","S2_score")

Baiting_df_final <- merge(Baiting_df,normalized_across_features,by = 0)
rownames(Baiting_df_final) <- Baiting_df_final$Row.names
Baiting_df_final$Row.names <- NULL

# These are the final LIBRA scores, adding these scores to the Seurat object
BC@meta.data <- merge(BC@meta.data,Baiting_df_final, by=0)
rownames(BC@meta.data) <- BC$Row.names
BC$Row.names <- NULL
remove(normalized_across_features,Baiting_df_final,Baiting_df,melted_Baiting_df)

# classify antigens
BC$antigen <- "Ag-"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Positive" & BC$cor_BF7_RBD.classification == "Negative" & BC$cor_S2.classification == "Negative"] <- "Ancestral RBD"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Negative" & BC$cor_BF7_RBD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "BF7 RBD"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Positive" & BC$cor_BF7_RBD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "Cross RBD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Positive" & BC$cor_BF7_NTD.classification == "Negative" & BC$cor_S2.classification == "Negative"] <- "Ancestral NTD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Negative" & BC$cor_BF7_NTD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "BF7 NTD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Positive" & BC$cor_BF7_NTD.classification == "Positive" & BC$cor_S2.classification == "Negative"] <- "Cross NTD"
BC$antigen[BC$cor_Ancestral_RBD.classification == "Negative" & BC$cor_BF7_RBD.classification == "Negative" & 
             BC$cor_Ancestral_NTD.classification == "Negative" & BC$cor_BF7_NTD.classification == "Negative" & BC$cor_S2.classification == "Positive"] <- "S2"

BC$antigen[BC$cor_Ancestral_RBD.classification == "Positive" & BC$cor_BF7_RBD.classification == "Negative" & BC$cor_S2.classification == "Positive" & BC$Ancestral_RBD_score < BC$S2_score] <- "S2"

BC$antigen[BC$cor_Ancestral_NTD.classification == "Positive" & BC$cor_BF7_NTD.classification == "Positive" & BC$cor_S2.classification == "Positive" & BC$Ancestral_NTD_score > BC$S2_score] <- "Cross NTD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Positive" & BC$cor_BF7_NTD.classification == "Negative" & BC$cor_S2.classification == "Positive" & BC$Ancestral_NTD_score > BC$S2_score] <- "Ancestral NTD"
BC$antigen[BC$cor_Ancestral_NTD.classification == "Negative" & BC$cor_BF7_NTD.classification == "Positive" & BC$cor_S2.classification == "Positive" & BC$BF7_NTD_score < BC$S2_score] <- "S2"

table(BC$bait.positive,BC$antigen)

#### 03.4 Batch2: add BCR info, classify isotype ####
BCR_batch2 <- read.table("./cellranger/batch2/vdj/filtered_contig.tsv",header = T,sep = "\t")
BCR_batch2$cell_id <- paste0(sapply(strsplit(BCR_batch2$sequence_id,"_"),'[',1),"_batch2")
BCR_anno <- read.csv("./cellranger/batch2/vdj/filtered_contig_annotations.csv")
BCR_anno <- BCR_anno[,c(3,27,28)]
colnames(BCR_anno) <- c("sequence_id","reads","umis")

BCR <- merge(BCR_batch2,BCR_anno,by = "sequence_id",all.x = T)
rm(BCR_batch2,BCR_anno)

BCR <- BCR %>%
  mutate(chain = ifelse(locus == 'IGH', 'Heavy', 'Light')) %>%
  group_by(cell_id) %>%
  mutate(
    type = case_when(
      sum(chain == "Heavy") >= 2 ~ "double",
      sum(chain == "Light") >= 2 ~ "double",
      sum(chain == "Heavy") == 0 & sum(chain == "Light") < 2 ~ "noHeavy",
      sum(chain == "Light") == 0 & sum(chain == "Heavy") < 2 ~ "noLight",
      TRUE ~ "single"
    )
  ) %>%
  ungroup()
table(BCR$type)

# double noHeavy noLight  single 
# 10469    6763     318   30826

BCR <- BCR %>%
  mutate(chain = ifelse(locus == 'IGH', 'Heavy', 'Light')) %>%
  group_by(cell_id, chain) %>%
  arrange(cell_id, desc(productive), desc(umis), desc(reads)) %>%
  slice_head(n=1)

BCR_H <- BCR[BCR$chain == "Heavy",]
BCR_H$SHM_count <- round((BCR_H$v_sequence_end - BCR_H$v_sequence_start +1) * (1 - (BCR_H$v_identity)/100))
BCR_H$SHM_freq <- 1 - BCR_H$v_identity/100
BCR_H$cdr3_length <- BCR_H$cdr3_end - BCR_H$cdr3_start + 1
colnames(BCR_H)[colnames(BCR_H) != "cell_id"] <- paste0("BCR_Heavy_", colnames(BCR_H)[colnames(BCR_H) != "cell_id"])

BCR_L <- BCR[BCR$chain == "Light",]
BCR_L$SHM_count <- round((BCR_L$v_sequence_end - BCR_L$v_sequence_start +1) * (1 - (BCR_L$v_identity)/100))
BCR_L$SHM_freq <- 1 - BCR_L$v_identity/100
BCR_L$cdr3_length <- BCR_L$cdr3_end - BCR_L$cdr3_start + 1
colnames(BCR_L)[colnames(BCR_L) != "cell_id"] <- paste0("BCR_Light_", colnames(BCR_L)[colnames(BCR_L) != "cell_id"])
BCR_merge_batch2 <- merge(BCR_H ,BCR_L, by='cell_id', all.x = T)

metadata <- BC@meta.data
metadata <- merge(metadata,BCR_merge_batch2,by = "cell_id",all.x = T)
rownames(metadata) <- metadata$cell_id
order <- colnames(BC)
metadata <- metadata[order,]

BC@meta.data <- metadata

BC <- add_BCR_info(BC, classify_col = "BCR_Heavy_c_call", new_col = "sub_Isotype")
table(BC$sub_Isotype)

BC@meta.data$Isotype <- "unknown"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype %in% c("IGHG1", "IGHG2","IGHG3","IGHG4","IGHG4A")] <- "IGHG"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype %in% c("IGHA1","IGHA2")] <- "IGHA"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype == "IGHD"] <- "IGHD"
BC@meta.data$Isotype[BC@meta.data$sub_Isotype == "IGHM"] <- "IGHM"

table(BC$Isotype)

rm(BCR,BCR_H,BCR_L,BCR_merge_batch2,metadata)

BC$RNA_snn_res.0.5 <- NULL
saveRDS(BC,"./for_publication/result/BC_batch2.rds")


#### 04. Merge 2 batch datasets, celltype annotation, add information and visualization ####
#### 04.1 BC WNN ####
# merge 2 batch dataset
BC_batch1 <- readRDS("./for_publication/result/BC_batch1.rds") #load seurat obj preprocessed in section2
BC_batch2 <- readRDS("./for_publication/result/BC_batch2.rds") #load seurat obj preprocessed in section3

BC_batch1$batch <- "batch1"
BC_batch2$batch <- "batch2"

BC <- merge(BC_batch1,BC_batch2)
rm(BC_batch1,BC_batch2)

BC@assays$Cor_Baiting <- NULL
# remove ZJJ because of cancer and 4 dose vaccination
BC <- subset(BC, sampleid != "SLE_ZJJ")

# remove batch effect using Harmony
DefaultAssay(BC) <- "RNA"
BC <- NormalizeData(BC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(BC,50)
BC <- RunHarmony(BC,reduction = "pca",group.by.vars = "sampleid",reduction.save = "harmony")
BC <- FindNeighbors(BC, reduction = "harmony", dims = 1:30)
BC <- FindClusters(BC, resolution = 0.5)
BC <- RunUMAP(BC, reduction = "harmony", dims = 1:30,reduction.name = "umap")
DimPlot(BC, reduction = 'umap', label = F, label.size = 6, pt.size = 0.8,group.by = c("sampleid"))+facet_wrap(~sampleid,ncol = 4)
DimPlot(BC, reduction = 'umap', label = T, label.size = 6, pt.size = 0.8,group.by = c("seurat_clusters"))
VlnPlot(BC, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), ncol = 2,pt.size = 0,group.by = "seurat_clusters")

My_DotPlot(BC,features = c(marker),group_by = "seurat_clusters",dot_scale = 6,coord_flip = FALSE)
My_DotPlot(BC,features = c("CD19","FCER2","IL4R","CR2","CD27","CD38","TFRC","CXCR5","ITGAX","FCRL5","SDC1","PRDM1"),group_by = "seurat_clusters",dot_scale = 6,coord_flip = FALSE)

# wnn, use ADT features(except IgG) for dimensional reduction
DefaultAssay(BC) <- 'Protein'
VariableFeatures(BC) <- c("CD21","FCRL5.1","CD71","CXCR5.1","IgM","CD11C","CD27.1","IgD","CD38.1")
BC <- NormalizeData(BC, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

BC <- FindMultiModalNeighbors(
  BC, reduction.list = list("harmony", "apca"),
  dims.list = list(1:30, 1:8), modality.weight.name = c("RNA.weight","Protein.weight"))
BC <- RunUMAP(BC, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
BC <- FindClusters(BC, graph.name = "wsnn", algorithm = 1, resolution = 1.8, verbose = FALSE)
DimPlot(BC, reduction = 'wnn.umap', label = T, label.size = 6, pt.size = 0.8,group.by = c("seurat_clusters"))+ggtitle("wnn_UMAP")+center.title()
DimPlot(BC, reduction = 'wnn.umap', label = F, label.size = 4, pt.size = 0.8,group.by = c("sampleid"))+facet_wrap(~sampleid,ncol = 3) + NoLegend()

DefaultAssay(BC) <- 'Protein'
FeaturePlot(BC,features = c("CD27.1","IgD","CD21","CD71","IgM","CD11C","FCRL5.1","CXCR5.1","CD38.1","IgG"),pt.size = 0.2,ncol = 5,cols = c("lightgrey","red4"))

My_DotPlot(BC,features = c("rna_FCER2","rna_IL4R","protein_CD27.1","rna_CD27","protein_IgD","protein_IgM","protein_IgG","protein_CD21","protein_CD38.1","protein_CD71","rna_TFRC","protein_CXCR5.1","protein_CD11C","rna_ITGAX","protein_FCRL5.1","rna_FCRL5","rna_PRDM1"),group_by = "seurat_clusters",dot_scale = 6,coord_flip = T)
VlnPlot(BC, features = c("rna_FCER2","rna_IL4R","protein_CD27.1","rna_CD27","protein_IgD","protein_IgM"), ncol = 2,pt.size = 0,group.by = "seurat_clusters")

# reannotation
BC$mainType <- "BC"
BC$subType <- ""
BC$subType[BC$seurat_clusters %in% c(1,5,6,7,9,10,11,12,13,15,19,20,21,23,24,27)] <- "rNAV"
BC$subType[BC$seurat_clusters %in% c(16)] <- "aNAV"
BC$subType[BC$seurat_clusters %in% c(16) & BC$Isotype %in% c("IGHG","IGHA")] <- "DN2"
BC$subType[BC$seurat_clusters %in% c(26)] <- "PB"
BC$subType[BC$seurat_clusters %in% c(17,18) & BC$Isotype != "IGHD"] <- "DN2"
BC$subType[BC$seurat_clusters %in% c(17,18) & BC$Isotype == "IGHD"] <- "aNAV"
BC$subType[BC$seurat_clusters %in% c(0,4) & BC$Isotype == "IGHD"] <- "rNAV"
BC$subType[BC$seurat_clusters %in% c(0,4) & BC$Isotype != "IGHD"] <- "uMBC"
BC$subType[BC$seurat_clusters %in% c(14)] <- "aMBC"
BC$subType[BC$seurat_clusters %in% c(25)] <- "uMBC"
BC$subType[BC$seurat_clusters %in% c(25) & BC$Isotype == "IGHD"] <- "rNAV"
BC$subType[BC$seurat_clusters %in% c(25) & BC$Isotype %in% c("IGHA","IGHG")] <- "sMBC"
BC$subType[BC$seurat_clusters %in% c(22)] <- "DN1"
BC$subType[BC$seurat_clusters %in% c(2,3,8)] <- "sMBC"
BC$subType[BC$seurat_clusters %in% c(2,3,8) & BC$Isotype %in% c("IGHM")] <- "uMBC"
BC$subType[BC$seurat_clusters %in% c(2,3,8) & BC$Isotype %in% c("IGHD")] <- "rNAV"
BC$subType[BC$subType == "uMBC" & BC$Isotype %in% c("IGHA","IGHG")] <- "sMBC"
BC$subType[BC$subType == "aMBC" & BC$Isotype %in% c("IGHD")] <- "rNAV"
BC$subType[BC$subType == "DN1" & BC$Isotype %in% c("IGHD")] <- "rNAV"
BC$subType[BC$subType == "rNAV" & BC$Isotype %in% c("IGHA","IGHG")] <- "sMBC"


BC$CD71 <- BC@assays$Protein@data["CD71",]
BC$CD27 <- BC@assays$Protein@data["CD27.1",]
BC$subType[BC$subType == "DN2" & BC$CD71 > 1.5] <- "aMBC" #we annotate CD71(hi) DN2 to aMBC
table(BC$subType,BC$Isotype)
BC$CD71 <- NULL
BC$CD27 <- NULL
BC$subType <- factor(BC$subType,levels = c("rNAV","aNAV","DN1","DN2","uMBC","sMBC","aMBC","PB"))

#### 04.2 Add metainfo of annotated B cells ####
BC$baiting <- "Ag+"
BC$baiting[BC$antigen == "Ag-"] <- "Ag-"
table(BC$baiting,BC$antigen)

BC$orig.ident <- "H"
BC$orig.ident[BC$sampleinfo %in% c("SLEI","SLEVI")] <- "SLE"

BC$sample_time <- paste0(BC$sampleinfo,"_",BC$timepoint)
table(BC$sample_time,BC$baiting)

BC$antigen_main <- "Ag-"
BC$antigen_main[BC$antigen %in% c("BF7 RBD","BF7 NTD")] <- "BF7.Specific"
BC$antigen_main[BC$antigen %in% c("Cross RBD","Cross NTD")] <- "Cross"
BC$antigen_main[BC$antigen %in% c("Ancestral RBD","Ancestral NTD")] <- "Ancestral"
BC$antigen_main[BC$antigen == "S2"] <- "S2"
table(BC$antigen_main,BC$antigen)

BC$epitope <- "Ag-"
BC$epitope[BC$antigen %in% c("BF7 RBD","Cross RBD","Ancestral RBD")] <- "RBD"
BC$epitope[BC$antigen %in% c("BF7 NTD","Cross NTD","Ancestral NTD")] <- "NTD"
BC$epitope[BC$antigen == "S2"] <- "S2"

saveRDS(BC,file = "./for_publication/result/BC_merged.rds")

#### 04.3 Visualization of annotated B cells ####

# Dimplot of WNN UMAP of total B cells
pdf("./figure/Fig4E.BC_WNN.pdf",width = 5.8,height = 5)
DimPlot(BC, reduction = 'wnn.umap', label = F, label.size = 6, pt.size = 0.8,group.by = c("subType"),cols = BC_col)+ggtitle("")+center.title() + NoLegend()
dev.off()

# Dotplot of B cell markers from protein level and RNA level
BC_sub <- subset(BC, subType != "PB")
BC_sub$subType <- factor(BC_sub$subType,levels = c("rNAV","aNAV","DN1","DN2","uMBC","sMBC","aMBC"))
pdf("./figure/Fig4F.BC_Dotplot_marker.pdf",width = 6.2,height = 4.8)
My_DotPlot(BC_sub,features = c("rna_FCER2","rna_IL4R","protein_CD27.1","protein_IgD",
                               "protein_IgM","protein_IgG","protein_CD21","protein_CD38.1",
                               "protein_CD71","protein_CXCR5.1","protein_CD11C","rna_ITGAX",
                               "protein_FCRL5.1","rna_FCRL5"),group_by = "subType",dot_scale = 6,coord_flip = T)
dev.off()

# VlnPlot of B cell markers of protein level and RNA level
DefaultAssay(BC) <- "Protein"
p <- VlnPlot(BC, features = c("IgD","CD27.1","IgM","CD21","CD71","FCRL5.1","CD11C","CXCR5.1","CD38.1","rna_PRDM1"), ncol = 5,pt.size = 0,group.by = "subType",cols = BC_col)
ggsave(p, filename = "./figure/FigS5B.Vlnplot_BC_protein_expr.pdf", width = 20, height = 5)

# Visualize Ag+ cells in different samples on wnnUMAP
Agneg <- BC@meta.data[BC$baiting == "Ag-","cell_id"]
Agpos_HI_early <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "HI_early","cell_id"]
Agpos_HVI_early <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "HVI_early","cell_id"]
Agpos_SLEI_early <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "SLEI_early","cell_id"]
Agpos_SLEVI_early <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "SLEVI_early","cell_id"]
Agpos_HI_late <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "HI_late","cell_id"]
Agpos_HVI_late <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "HVI_late","cell_id"]
Agpos_SLEI_late <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "SLEI_late","cell_id"]
Agpos_SLEVI_late <- BC@meta.data[BC$baiting == "Ag+" & BC$sample_time == "SLEVI_late","cell_id"]

metadata <- BC@meta.data[,c("cell_id","antigen")]

cell_list <- list(Agpos_HI_early,Agpos_HVI_early,Agpos_SLEI_early,Agpos_SLEVI_early,Agpos_HI_late,Agpos_HVI_late,Agpos_SLEI_late,Agpos_SLEVI_late)
names(cell_list) <- c("HI_early","HVI_early","SLEI_early","SLEVI_early","HI_late","HVI_late","SLEI_late","SLEVI_late")

p_main <- list()

for (i in 1:length(cell_list)) {
  sample <- subset(BC, cell_id %in% c(Agneg, cell_list[[i]]))
  umap_coords <- as.data.frame(Embeddings(sample, reduction = "wnn.umap"))
  colnames(umap_coords) <- c("wnnUMAP_1", "wnnUMAP_2")
  umap_coords$cell_id <- rownames(umap_coords)
  umap_coords <- merge(umap_coords, metadata, by = "cell_id")
  
  umap_agneg <- umap_coords %>% filter(antigen == "Ag-")
  umap_agpos <- umap_coords %>% filter(antigen != "Ag-")
  
  x_range <- range(umap_coords$wnnUMAP_1, na.rm = TRUE)
  y_range <- range(umap_coords$wnnUMAP_2, na.rm = TRUE)
  
  main_plot <- ggplot() +
    geom_point(data = umap_agneg, aes(wnnUMAP_1, wnnUMAP_2),
               color = "grey95", size = 0.2) +
    geom_point(data = umap_agpos, aes(wnnUMAP_1, wnnUMAP_2, color = antigen),
               size = 2.2) +
    scale_color_manual(values = Antigen) +
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    scale_y_continuous(limits = y_range, expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black", size = 1),
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    )
  
  top_hist <- ggplot(umap_agpos, aes(wnnUMAP_1)) +
    geom_histogram(fill = "#EEC591", color = "black", alpha = 0.5, 
                   bins = 50, boundary = x_range[1], closed = "left",
                   linewidth = 0.8) +
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
  
  right_hist <- ggplot(umap_agpos, aes(wnnUMAP_2)) +
    geom_histogram(fill = "#EEC591", color = "black", alpha = 0.5, 
                   bins = 50, boundary = y_range[1], closed = "left",
                   linewidth = 0.8) +
    scale_x_continuous(limits = y_range, expand = c(0, 0)) +
    coord_flip() +
    theme_void() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
  
  p_main[[i]] <- plot_grid(
    top_hist, NULL,
    main_plot, right_hist,
    ncol = 2, nrow = 2,
    rel_widths = c(4, 1),
    rel_heights = c(1, 4),
    align = "hv",
    axis = "tblr"
  )
}

p_main_with_margin <- lapply(p_main, function(p) {
  ggdraw(p) + theme(plot.margin = margin(5, 5, 5, 5, "mm"))
})

pdf("Fig3D.AgBC_WNN_UMAP.pdf", width = 22, height = 9)
grid.arrange(grobs = p_main_with_margin,
             ncol = 4,
             nrow = 2)
dev.off()

rm(Agneg,Agpos_HI_early,Agpos_HVI_early,Agpos_SLEI_early,Agpos_SLEVI_early,Agpos_HI_late,
   Agpos_HVI_late,Agpos_SLEI_late,Agpos_SLEVI_late,cell_list,p_main,p_main_with_margin,
   sample,umap_coords,umap_agneg,umap_agpos,top_hist,right_hist,main_plot,x_range,y_range)

# Distribution of LIBRA-seq score of Ag+ B cells
BC$Ancestral_RBD_score[is.na(BC$Ancestral_RBD_score)] <- 0
BC$Ancestral_NTD_score[is.na(BC$Ancestral_NTD_score)] <- 0
BC$BF7_RBD_score[is.na(BC$BF7_RBD_score)] <- 0
BC$BF7_NTD_score[is.na(BC$BF7_NTD_score)] <- 0
BC$S2_score[is.na(BC$S2_score)] <- 0

p <- list()
p[[1]] <- FeatureScatter(BC, feature1 = "Ancestral_RBD_score", feature2 = "BF7_RBD_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("Ancestral RBD scores") + ylab("BF.7 RBD scores")+NoLegend()
p[[2]] <- FeatureScatter(BC, feature1 = "Ancestral_NTD_score", feature2 = "BF7_NTD_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("Ancestral NTD scores") + ylab("BF.7 NTD scores")+NoLegend()
p[[3]] <- FeatureScatter(BC, feature1 = "BF7_RBD_score", feature2 = "S2_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("BF.7 RBD scores") + ylab("S2 scores")+NoLegend()
p[[4]] <- FeatureScatter(BC, feature1 = "BF7_NTD_score", feature2 = "S2_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("BF.7 NTD scores") + ylab("S2 scores")+NoLegend()
p[[5]] <- FeatureScatter(BC, feature1 = "Ancestral_RBD_score", feature2 = "S2_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("Ancestral RBD scores") + ylab("S2 scores")+NoLegend()
p[[6]] <- FeatureScatter(BC, feature1 = "Ancestral_NTD_score", feature2 = "S2_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("Ancestral NTD scores") + ylab("S2 scores")+NoLegend()

p[[7]] <- FeatureScatter(BC, feature1 = "BF7_RBD_score", feature2 = "BF7_NTD_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("BF7 RBD scores") + ylab("BF7 NTD scores")+NoLegend()
p[[8]] <- FeatureScatter(BC, feature1 = "BF7_RBD_score", feature2 = "Ancestral_NTD_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("BF7 RBD scores") + ylab("Ancestral NTD scores")+NoLegend()
p[[9]] <- FeatureScatter(BC, feature1 = "Ancestral_RBD_score", feature2 = "BF7_NTD_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("Ancestral RBD scores") + ylab("BF7 NTD scores")+NoLegend()
p[[10]] <- FeatureScatter(BC, feature1 = "Ancestral_RBD_score", feature2 = "Ancestral_NTD_score", pt.size = 1, group.by = "antigen", cols = Antigen) +
  ggtitle("") + xlab("Ancestral RBD scores") + ylab("Ancestral NTD scores")+NoLegend()

pdf("FigS4D.Ag+_LIBRA_score.pdf",width = 10,height = 6)
grid.arrange(grobs = p, ncol = 3, nrow = 2)
dev.off()


# Ag+ MBC cells epitope distribution
metadata <- BC@meta.data
metadata <- metadata[metadata$baiting == "Ag+" & metadata$subType %in% c("DN1","uMBC","sMBC","DN2","aMBC"),]
table(metadata$sample_time)

tt = as.data.frame(round(100*prop.table(x = table(metadata$epitope,metadata$sample_time),margin = 2),2))
names(tt) = c("antigen","sample", "Freq")
tt$antigen <- factor(tt$antigen,levels = c("RBD","NTD","S2"))
tt$sample <- factor(tt$sample,levels = c("HI_early","HI_late","HVI_early","HVI_late","SLEI_early","SLEI_late","SLEVI_early","SLEVI_late"))

p <- ggplot(tt,aes(y = Freq, x = factor(sample), fill = antigen)) + geom_bar(stat="identity")+
  scale_fill_manual(values = c("#FA8090", "#f9ab60", "#A6D0E6")) + 
  labs(x="", y="Spike subdomain+ MBC (%)") +
  theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,colour="black", size = rel(2)),axis.line.x = element_line(colour = "black", size = 1),
        axis.text.y=element_text(colour="black", size = rel(2)),legend.position = "None",axis.line.y = element_line(colour = "black", size = 1),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),legend.text = element_text(size = 15),panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA))

ggsave("FigS5A.Ag+_MBC_distribution_byepitope.pdf",p,width = 6,height = 5.5)


# Ag+MBCs epitope distribution

metadata <- BC@meta.data[BC$antigen_main %in% c("BF7.Specific","Cross") & BC$subType %in% c("DN1","uMBC","sMBC","DN2","aMBC"),] 

table(metadata$sample_time)

tt = as.data.frame(round(100*prop.table(x = table(metadata$antigen_main,metadata$sample_time),margin = 2),2))
names(tt) = c("antigen","sample", "Freq")
tt$antigen <- factor(tt$antigen,levels = c("BF7.Specific","Cross"))
tt$sample <- factor(tt$sample,levels = c("HI_early","HI_late","HVI_early","HVI_late","SLEI_early","SLEI_late","SLEVI_early","SLEVI_late"))
p <- ggplot(tt,aes(y = Freq, x = factor(sample), fill = antigen)) + geom_bar(stat="identity")+
  scale_fill_manual(values = c("#88C0C0","grey90")) +
  labs(x="", y="Spike subdomain+ MBC (%)") +
  theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,colour="black", size = rel(2)),axis.line.x = element_line(colour = "black", size = 1),
        axis.text.y=element_text(colour="black", size = rel(2)),legend.position = "None",axis.line.y = element_line(colour = "black", size = 1),
        axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),legend.text = element_text(size = 15),panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA))

ggsave("FigS5B.Imprinting_longitudinal.pdf",p,width = 6,height = 5.5)

rm(metadata,tt)

#### 05. Transcriptome downstream analysis ####
# When calculate DEG and pathway enrichment, I added Ig genes back to RNA expression matrix
raw <- readRDS("raw_combined.rds")
common_cells <- intersect(colnames(BC@assays$RNA), colnames(BC@assays$RNA))
BC[["RNA"]] <- CreateAssayObject(raw@assays$RNA[, common_cells])
rm(raw,common_cells)

DefaultAssay(BC) <- "RNA"
BC <- NormalizeData(BC)

#### 05.1 DEG analysis and visualization ####
#heatmap to visualize DEGs
gene_list <- list(
  "BCR signaling and B cell activation" = c("CD79B", "MS4A1", "CD79A", "CD27", "SYK"),
  "Antigen processing and presentation" = c("HLA-DRB1", "CD74","HLA-DOB"),
  "Inflammatory response" = c("CD69", "NFKBIA", "TXNIP", "PELI1", "KLF6", "IFI44L", "ISG20"),
  "Glucocorticoid stimulus" = c("TSC22D3", "DUSP1")
)
genes_to_plot <- unlist(gene_list)
gene_groups <- data.frame(
  gene = genes_to_plot,
  pathway = rep(names(gene_list), sapply(gene_list, length))
)
rownames(gene_groups) <- gene_groups$gene

BC_heat <- ScaleData(BC, features = genes_to_plot)
DefaultAssay(BC_heat) <- "RNA"
AgBC <- subset(BC_heat,sampleinfo %in% c("SLEVI","HVI") & subType %in% c("DN2","DN1","uMBC","sMBC","aMBC") & baiting == "Ag+")

expr_matrix <- GetAssayData(AgBC, slot = "scale.data")[genes_to_plot, ]

cell_metadata <- AgBC@meta.data[, "sampleinfo", drop = FALSE]
cell_order <- order(cell_metadata$sampleinfo)
expr_matrix_ordered <- expr_matrix[, cell_order]

sample_colors <- setNames(
  c(HVI,SLEVI), 
  unique(cell_metadata$sampleinfo)
)

col_annotation <- HeatmapAnnotation(
  Sample = cell_metadata$sampleinfo[cell_order],
  col = list(Sample = sample_colors),
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10)
)

row_split <- factor(gene_groups[genes_to_plot, "pathway"], 
                    levels = names(gene_list))

col_fun <- colorRamp2(c(-2, 0, 2), c("#0571b0", "white", "#ca0020"))

row_annotation <- rowAnnotation(
  Pathway = gene_groups[genes_to_plot, "pathway"],
  col = list(Pathway = setNames(
    c("#5242C3", "#074070","#D73114","#FF8C00"),
    names(gene_list)
  )),
  show_legend = FALSE,
  width = unit(5, "mm")
)

pdf("Fig3E.Heatmap_DEG_grouped.pdf", width = 9, height = 4)
ht <- Heatmap(
  expr_matrix_ordered,
  name = "Expression",
  col = col_fun,
  
  top_annotation = col_annotation,
  #right_annotation = row_annotation,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = cell_metadata$sampleinfo[cell_order],
  column_title_gp = gpar(fontsize = 12),
  column_title_rot = 0,
  column_gap = unit(2, "mm"),
  
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  row_split = row_split,
  row_title_gp = gpar(fontsize = 11),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),
  
  heatmap_legend_param = list(
    title = "Expression",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontsize = 10)
  ),border = TRUE
)

draw(ht, heatmap_legend_side = "right")
dev.off()


#heatmap to visualize DEGs
gene_list <- list(
  "BCR signaling and B cell activation" = c("CD79B", "MS4A1", "CD79A", "SYK"),
  "Antigen processing and presentation" = c("HLA-DRB1", "CD74","HLA-DOB"),
  "Inflammatory response" = c("CD69", "NFKBIA", "TXNIP", "PELI1", "KLF6", "IFI44L", "ISG20"),
  "Glucocorticoid stimulus" = c("TSC22D3", "DUSP1")
)
genes_to_plot <- unlist(gene_list)
gene_groups <- data.frame(
  gene = genes_to_plot,
  pathway = rep(names(gene_list), sapply(gene_list, length))
)
rownames(gene_groups) <- gene_groups$gene

BC_heat <- ScaleData(BC, features = genes_to_plot)
DefaultAssay(BC_heat) <- "RNA"
NAV <- subset(BC_heat,sampleinfo %in% c("SLEVI","HVI") & subType %in% c("rNAV","aNAV"))

expr_matrix <- GetAssayData(NAV, slot = "scale.data")[genes_to_plot, ]

cell_metadata <- NAV@meta.data[, "sampleinfo", drop = FALSE]
group_order <- c("HVI", "SLEVI")

set.seed(42)
idx_list <- lapply(group_order, function(g) {
  sample(which(cell_metadata$sampleinfo == g))
})
cell_order <- unlist(idx_list)

expr_matrix_ordered <- expr_matrix[, cell_order]

sample_colors <- setNames(
  c(HVI,SLEVI), 
  unique(cell_metadata$sampleinfo)
)

col_annotation <- HeatmapAnnotation(
  Sample = cell_metadata$sampleinfo[cell_order],
  col = list(Sample = sample_colors),
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10)
)

row_split <- factor(gene_groups[genes_to_plot, "pathway"], 
                    levels = names(gene_list))

col_fun <- colorRamp2(c(-2, 0, 2), c("#0571b0", "white", "#ca0020"))

row_annotation <- rowAnnotation(
  Pathway = gene_groups[genes_to_plot, "pathway"],
  col = list(Pathway = setNames(
    c("#5242C3", "#074070","#D73114","#FF8C00"),
    names(gene_list)
  )),
  show_legend = FALSE,
  width = unit(5, "mm")
)

pdf("FigS5C.Heatmap_NAV_DEG_grouped.pdf", width = 9, height = 4)
ht <- Heatmap(
  expr_matrix_ordered,
  name = "Expression",
  col = col_fun,
  
  top_annotation = col_annotation,
  #right_annotation = row_annotation,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = cell_metadata$sampleinfo[cell_order],
  column_title_gp = gpar(fontsize = 12),
  column_title_rot = 0,
  column_gap = unit(2, "mm"),
  
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  row_split = row_split,
  row_title_gp = gpar(fontsize = 11),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),
  
  heatmap_legend_param = list(
    title = "Expression",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontsize = 10)
  ),border = TRUE
)

draw(ht, heatmap_legend_side = "right")
dev.off()


rm(BC_heat,row_annotation,col_fun,row_split,col_annotation,
   sample_colors,cell_metadata,cell_order,expr_matrix_ordered,
   expr_matrix,AgBC,gene_groups,genes_to_plot,ht,gene_list)

#### 05.2 Add module score and GSEA ####
# Loading the gene sets from the MSigDB database
m_df<- msigdbr(species = "Homo sapiens", category = "H") #hallmark gene sets
m_df2<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG") #curated gene sets
m_df3<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME") #curated gene sets
m_df4<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP") #ontology gene sets

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets2 <- m_df2 %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets3 <- m_df3 %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets4 <- m_df4 %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets <- c(fgsea_sets,fgsea_sets2,fgsea_sets3,fgsea_sets4)
rm(m_df,m_df2,m_df3,m_df4,fgsea_sets2,fgsea_sets3,fgsea_sets4)

DEG <- FindMarkers(AgBC, ident.1 = "SLE",ident.2 = "H",group.by = "orig.ident",logfc.threshold = 0.1)
DEG$gene <- rownames(DEG)
#DEG_Ag_MBC <- DEG[DEG$p_val_adj < 0.05,]
write.xlsx(DEG,file = "DEG_Ag+MBC_early_late_SLEVIvsHVI.xlsx")

preranked_list <- DEG %>% 
  mutate(ranking_metric = avg_log2FC) %>%
  arrange(desc(ranking_metric)) %>% 
  dplyr::select(gene, ranking_metric) %>% deframe

fgseaResults <- fgsea(pathways = fgsea_sets, stats = preranked_list) %>%
  as_tibble() %>% arrange(desc(NES))
fgseaResults_filter <- fgseaResults %>% filter(padj < 0.05)
write.xlsx(fgseaResults_filter,file = "GSEA_Ag+MBC_early_late_SLEVIvsHVI.xlsx")

p <- plotEnrichment(fgsea_sets[["HALLMARK_INFLAMMATORY_RESPONSE"]],
                    preranked_list) + labs(title="Inflammatory response")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("Fig3F.GSEA_HALLMARK_INFLAMMATORY_RESPONSE.pdf",p,width = 4,height = 2.5)

p <- plotEnrichment(fgsea_sets[["GOBP_B_CELL_ACTIVATION"]],
                    preranked_list) + labs(title="B cell activation")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("Fig3F.GSEA_GOBP_B_CELL_ACTIVATION.pdf",p,width = 4,height = 2.5)

p <- plotEnrichment(fgsea_sets[["KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"]],
                    preranked_list) + labs(title="Antigen processing and presentation")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("Fig3F.GSEA_KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION.pdf",p,width = 4,height = 2.5)

p <- plotEnrichment(fgsea_sets[["REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS"]],
                    preranked_list) + labs(title="BCR signaling")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("Fig3F.GSEA_REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR.pdf",p,width = 4,height = 2.5)

#analyze total Naive B cell
DEG <- FindMarkers(NAV, ident.1 = "SLE",ident.2 = "H",group.by = "orig.ident",logfc.threshold = 0.1)
DEG$gene <- rownames(DEG)
DEG_Ag_NAV <- DEG[DEG$p_val_adj < 0.05,]
write.xlsx(DEG,file = "DEG_NAV_early_late_SLEVIvsHVI.xlsx")

preranked_list <- DEG %>% 
  mutate(ranking_metric = avg_log2FC) %>%
  arrange(desc(ranking_metric)) %>% 
  dplyr::select(gene, ranking_metric) %>% deframe

fgseaResults <- fgsea(pathways = fgsea_sets, stats = preranked_list) %>%
  as_tibble() %>% arrange(desc(NES))
fgseaResults_filter <- fgseaResults %>% filter(padj < 0.05)
write.xlsx(fgseaResults_filter,file = paste0("GSEA_Ag+NAV_early_late_SLEVIvsHVI.xlsx"))

p <- plotEnrichment(fgsea_sets[["HALLMARK_INFLAMMATORY_RESPONSE"]],
                    preranked_list) + labs(title="Inflammatory response")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("FigS5D.GSEA_NAV_HALLMARK_INFLAMMATORY_RESPONSE.pdf",p,width = 4,height = 2.5)

p <- plotEnrichment(fgsea_sets[["GOBP_B_CELL_ACTIVATION"]],
                    preranked_list) + labs(title="B cell activation")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("FigS5D.GSEA_NAV_GOBP_B_CELL_ACTIVATION.pdf",p,width = 4,height = 2.5)

p <- plotEnrichment(fgsea_sets[["REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS"]],
                    preranked_list) + labs(title="BCR signaling")+labs(x="Rank",y="Enrichment score")+theme_cowplot()+center.title() + 
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
ggsave("FigS5D.GSEA_NAV_REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR.pdf",p,width = 4,height = 2.5)

#### 06. BCR repertoire ####
#### 06.1 Define clone ####
table(BC$sampleid)

BC$sample <- BC$sampleid
BC$sample[BC$sampleid %in% c("H002-1","H002-2")] <- "H002"
BC$sample[BC$sampleid %in% c("H033-1","H033-2")] <- "H033"
BC$sample[BC$sampleid %in% c("H050-1","H050-2")] <- "H050"
BC$sample[BC$sampleid %in% c("H052-1","H052-2")] <- "H052"
BC$sample[BC$sampleid %in% c("R014-1")] <- "R014"
BC$sample[BC$sampleid %in% c("R017-1","R017-2")] <- "R017"
BC$sample[BC$sampleid %in% c("CR-3")] <- "CR-3"
BC$sample[BC$sampleid %in% c("CR-11")] <- "CR-11"
BC$sample[BC$sampleid %in% c("SLE_WL")] <- "SLE_WL"

table(BC$sample)

BCR <- BC@meta.data[is.na(BC$BCR_Heavy_sequence_id) == F,c(1,49:284)]


names(BCR) <- gsub("BCR_Heavy_", "", names(BCR))
write.table(BCR, './result/BCR/changeo/BCR_H.tsv', sep = '\t', quote = F, row.names = F)


library(shazam)
dist_ham <- distToNearest(BC@meta.data[!is.na(BC$BCR_Heavy_locus),], 
                          locusColumn="BCR_Heavy_locus",
                          sequenceColumn="BCR_Heavy_junction", 
                          vCallColumn="BCR_Heavy_v_call", jCallColumn="BCR_Heavy_j_call",
                          model="ham", normalize="len", nproc=1, fields="sample")

ggplot(subset(dist_ham, !is.na(dist_nearest)), 
       aes(x=dist_nearest)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.241, color="firebrick", linetype=2) +
  labs(x = "Grouped Hamming distance of Heavy chain", y = "Count") + 
  facet_grid(sample ~ ., scales="free_y") + 
  theme_bw()

output <- findThreshold(dist_ham$dist_nearest, method="gmm", model="gamma-gamma")
# Plot distance histogram, Gaussian fits, and optimum threshold
plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
print(output) # FIND THRESHOLD FOR 3 TIMES: 0.2409758, 0.2458953, 0.2414974, CHOOSE 0.241 AS THRESHOLD

dist_ham <- distToNearest(BC@meta.data[!is.na(BC$BCR_Light_locus),], 
                          locusColumn="BCR_Light_locus",locusValues = c("IGK","IGL"),
                          sequenceColumn="BCR_Light_junction", 
                          vCallColumn="BCR_Light_v_call", jCallColumn="BCR_Light_j_call",
                          model="ham", normalize="len", nproc=1, fields="sample")

ggplot(subset(dist_ham, !is.na(dist_nearest)), 
       aes(x=dist_nearest)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.07, color="firebrick", linetype=2) +
  labs(x = "Grouped Hamming distance of Light chain", y = "Count") + 
  facet_grid(sample ~ ., scales="free_y") + 
  theme_bw()


#define clone by changeo
#DefineClones.py -d BCR_H.tsv --act set --model ham --norm len --failed --format airr --dist 0.241 --outdir dir


metadata <- BC@meta.data

# We used heavy chain to define clone
clone_pass <- read.table("./result/BCR/changeo/BCR_H_clone-pass.tsv",header = T, sep = "\t")
clone_pass <- clone_pass[clone_pass$productive == "TRUE",]
rownames(clone_pass) <- clone_pass$cell_id
H_clone <- clone_pass[,c("cell_id","clone_id")]
metadata <- merge(metadata,H_clone,by = "cell_id",all.x = T)
metadata$clone_id[is.na(metadata$clone_id)] <- "unknown"
rownames(metadata) <- metadata$cell_id
order <- BC$cell_id
metadata <- metadata[order,]
BC@meta.data <- metadata

rm(clone_pass,H_clone,metadata,dist_ham,output,order,BCR)
saveRDS(BC,"BC_merged.rds")

#### 06.2 Clonal relationship ####

### for Ag+ IgG MBC and total Naive and IgM MBC relationship
samplename <- c("H002","H033","H050","H052","R014","R017","SLE_WL","CR-11","CR-3")
metadata <- BC@meta.data[BC$clone_id != "unknown",]
meta_NAV <- metadata[metadata$subType %in% c("rNAV","aNAV"),]
meta_IgM <- metadata[metadata$subType %in% c("DN1","DN2","uMBC","sMBC","aMBC") & metadata$Isotype == "IGHM",]
meta_IgG <- metadata[metadata$subType %in% c("DN1","DN2","uMBC","sMBC","aMBC") & metadata$Isotype == "IGHG" & metadata$baiting == "Ag+",]

metadata <- rbind(meta_NAV,meta_IgM,meta_IgG)
metadata$subType <- as.character(metadata$subType)
metadata$clone_id <- paste0(metadata$clone_id,"_",metadata$sample)
table(metadata$sampleinfo)

# prepare data for circos plot
title <- samplename[i] #i=1:9
clone_df <- metadata[metadata$sample == title,c("sample","timepoint","baiting","BCR_Heavy_SHM_count","Isotype","antigen","subType","clone_id")]
#use Circos_prepare.R to generate circos files for each individual

# generate circos files for H002
H002_cells <- read.csv("./H002/H002_cells.csv")

#select cells with relationships between different subTypes and visualize in circos plot
H002_cells_filtered <- H002_cells[,c(1:8)] %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1) %>%
  ungroup()
clone_df <- H002_cells_filtered

H002_links_filtered <- H002_cells %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1, any(antigen != "Ag-", na.rm = TRUE)) %>%
  ungroup()

H002_cells <- read.csv("./H002/H002_cells.csv")
timepoint <- H002_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "H002_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
Ag <- H002_cells[,c(9,10,11,3)]
Ag$baiting[Ag$baiting == "Ag+"] <- 1
Ag$baiting[Ag$baiting == "Ag-"] <- 0
write.table(Ag,file = "H002_Ag.txt",quote = F,row.names = F,col.names = F,sep = "\t")

rm(H002_cells,H002_link,timepoint)

# generate circos files for H050
H033_cells <- read.csv("./H033/H033_cells.csv")

H033_cells_filtered <- H033_cells[,c(1:8)] %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1) %>%
  ungroup()
clone_df <- H033_cells_filtered

H033_links_filtered <- H033_cells %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1, any(antigen != "Ag-", na.rm = TRUE)) %>%
  ungroup()

H033_cells <- read.csv("./H033/H033_cells.csv")
timepoint <- H033_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "H033_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
Ag <- H033_cells[,c(9,10,11,3)]
Ag$baiting[Ag$baiting == "Ag+"] <- 1
Ag$baiting[Ag$baiting == "Ag-"] <- 0
write.table(Ag,file = "H033_Ag.txt",quote = F,row.names = F,col.names = F,sep = "\t")

rm(H033_cells,timepoint,H033_cells_filtered,Ag)


# generate circos files for H050
H050_cells <- read.csv("./H050/H050_cells.csv")

H050_cells_filtered <- H050_cells[,c(1:8)] %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1) %>%
  ungroup()
clone_df <- H050_cells_filtered

H050_links_filtered <- H050_cells %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1, any(antigen != "Ag-", na.rm = TRUE)) %>%
  ungroup()

H050_cells <- read.csv("./H050/H050_cells.csv")
timepoint <- H050_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "H050_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
Ag <- H050_cells[,c(9,10,11,3)]
Ag$baiting[Ag$baiting == "Ag+"] <- 1
Ag$baiting[Ag$baiting == "Ag-"] <- 0
write.table(Ag,file = "H050_Ag.txt",quote = F,row.names = F,col.names = F,sep = "\t")

rm(H050_cells,Ag,timepoint,H050_cells_filtered)

# generate circos files for H052
H052_cells <- read.csv("./H052/H052_cells.csv")

H052_cells_filtered <- H052_cells[,c(1:8)] %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1) %>%
  ungroup()
clone_df <- H052_cells_filtered

H052_links_filtered <- H052_cells %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1, any(antigen != "Ag-", na.rm = TRUE)) %>%
  ungroup()

H052_cells <- read.csv("./H052/H052_cells.csv")
timepoint <- H052_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "H052_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
Ag <- H052_cells[,c(9,10,11,3)]
Ag$baiting[Ag$baiting == "Ag+"] <- 1
Ag$baiting[Ag$baiting == "Ag-"] <- 0
write.table(Ag,file = "H052_Ag.txt",quote = F,row.names = F,col.names = F,sep = "\t")

rm(H052_cells,Ag,H052_cells_filtered)

### for Total B cell relationship
samplename <- c("HI","HVI","SLEI","SLEVI")
metadata <- BC@meta.data[BC$clone_id != "unknown",]
metadata$subType <- as.character(metadata$subType)
metadata$clone_id <- paste0(metadata$clone_id,"_",metadata$sample)
table(metadata$sampleinfo)

# prepare data for circos plot

title <- samplename[i] #i=1:4

clone_df <- metadata[metadata$sampleinfo == title,c("sample","timepoint","baiting","BCR_Heavy_SHM_count","Isotype","antigen","subType","clone_id")]

# generate circos files for HI
HI_cells <- read.csv("./HI/HI_cells.csv")

SLEI_cells_filtered <- SLEI_cells[,c(1:8)] %>%
  group_by(clone_id) %>%
  filter(n() >= 2, n_distinct(subType) > 1) %>%
  ungroup()

timepoint <- HI_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "./HI/HI_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
rm(HI_cells,HI_link,timepoint)

# generate circos files for HVI
HVI_cells <- read.csv("./HVI/HVI_cells.csv")

timepoint <- HVI_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "./HVI/HVI_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
rm(HVI_cells,HVI_link,timepoint)

# generate circos files for SLEI
SLEI_cells <- read.csv("./SLEI/SLEI_cells.csv")

timepoint <- SLEI_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "./SLEI/SLEI_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
rm(SLEI_cells,SLEI_link,timepoint)

# generate circos files for SLEVI
SLEVI_link <- read.csv("./SLEVI/SLEVI_link.csv")
SLEVI_cells <- read.csv("./SLEVI/SLEVI_cells.csv")

timepoint <- SLEVI_cells[,c(9,10,11,2)]
timepoint$timepoint[timepoint$timepoint == "early"] <- 1
timepoint$timepoint[timepoint$timepoint == "late"] <- 0
write.table(timepoint,file = "./SLEVI/SLEVI_timepoint.txt",quote = F,row.names = F,col.names = F,sep = "\t")
rm(SLEVI_cells,SLEVI_link,timepoint)


#### 06.3 Lineage trees ####
library(dowser)
BC_sub <- subset(BC,cell_id %in% metadata$cell_id)

# Prepare BCR data for H002
clone_info <- BC_sub@meta.data[BC_sub$clone_id %in% c(11837) & BC_sub$sample == "H002",]
clone_info$subType <- as.character(clone_info$subType)
clone_info$antigen_main <- as.character(clone_info$antigen_main)
names(clone_info) <- gsub("BCR_Heavy_", "", names(clone_info))

clone <- unique(clone_info$clone_id)
for (i in 1:length(clone)){
  print(clone[i])
  sub_clone_info <- clone_info[clone_info$clone_id == clone[i],]
  for (j in 1:nrow(sub_clone_info)){
    title <- sub_clone_info[j,"subType"]
    seq <- sub_clone_info[j,"germline_alignment"]
    if (is.character(seq) && length(seq) == 1) {
      seq <- strsplit(seq, "")[[1]]
    }
    write.fasta(sequences = list(seq), names = title, file.out = paste0("./tree/H002_",clone[i],"_germ.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}

# assign germline alignment with cdr3(ACTGT...TGGGG) masked according to igblast

clone_info$germline_alignment_d_mask[clone_info$clone_id == 11837] <- "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTACGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG"
clones = formatClones(clone_info, traits=c("Isotype","subType","timepoint","antigen_main"))
print(clones)

save(clones,file = "./tree/H002_clone_info.RData")

# use linux R for igphyml
load("H002_clone_info.RData")
clones = getTrees(clones, nproc=1,build="igphyml",exec="./software/igphyml/src/igphyml")
save(clones,file = "H002_clone_tree.RData")

load("./tree/H002_clone_tree.RData")
clones = getTrees(clones, nproc=1)

custom_palette=c(BC_col,
                 "Germline"= "black",
                 "early" = "yellow",
                 "late" = "darkred",
                 "BF7.specific" = BF7,
                 "Cross" = Cross,
                 "Ancestral" = Ancestral,
                 "S2" = "#CDC9A5",
                 "Ag-" = "grey",
                 "IGHD" = "#A9A9A9",
                 "IGHA" = "#44AA99",
                 "IGHG" = "#ff7f00",
                 "IGHM" = "#BCD3DF",
                 "unknown" = "white"
)

plots = plotTrees(clones, tips = "subType",tipsize = 6,palette=custom_palette)
plots[[1]]
treesToPDF(plots, file="./tree/Fig6D.H002_trees_subType.pdf", nrow=1, ncol=1,width = 4,height=3)

plots_time = plotTrees(clones, tips = "timepoint",tipsize = 6,palette=custom_palette)
plots_time[[1]]
treesToPDF(plots_time, file="./tree/Fig6D.H002_trees_timepoint.pdf", nrow=1, ncol=1,width = 4,height=3)

plots_ag = plotTrees(clones, tips = "antigen_main",tipsize = 6,palette=custom_palette)
plots_ag[[1]]
treesToPDF(plots_ag, file="./tree/Fig6D.H002_trees_antigen.pdf", nrow=1, ncol=1,width = 4,height=3)

plots_Iso = plotTrees(clones, tips = "Isotype",tipsize = 6,palette=custom_palette)
plots_Iso[[1]]
treesToPDF(plots_Iso, file="./tree/Fig6D.H002_trees_Isotype.pdf", nrow=1, ncol=1,width = 4,height=3)

# Prepare BCR data for H050
clone_info <- BC_sub@meta.data[BC_sub$clone_id %in% c(721) & BC_sub$sample == "H050",]
clone_info$subType <- as.character(clone_info$subType)
clone_info$antigen_main <- as.character(clone_info$antigen_main)
names(clone_info) <- gsub("BCR_Heavy_", "", names(clone_info))

clone <- unique(clone_info$clone_id)
for (i in 1:length(clone)){
  print(clone[i])
  sub_clone_info <- clone_info[clone_info$clone_id == clone[i],]
  for (j in 1:nrow(sub_clone_info)){
    title <- sub_clone_info[j,"subType"]
    seq <- sub_clone_info[j,"germline_alignment"]
    if (is.character(seq) && length(seq) == 1) {
      seq <- strsplit(seq, "")[[1]]
    }
    write.fasta(sequences = list(seq), names = title, file.out = paste0("./for_publication/figure/tree/H050_",clone[i],"_germ.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}

clone_info$germline_alignment_d_mask[clone_info$clone_id == 721] <- "CAGATCACCTTGAAGGAGTCTGGTCCTACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCACCTTCTCTGGGTTCTCACTCAGCACTAGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGGCCCTGGAGTGGCTTGCACTCATTTATTGGGATGATGATAAGCGCTACAGCCCATCTCTGAAGAGCAGGCTCACCATCACCAAGGACACCTCCAAAAACCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGACACAGCCACATATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clones = formatClones(clone_info, traits=c("Isotype","subType","timepoint","antigen_main"))

print(clones)

save(clones,file = "./tree/H050_clone_info.RData")

# use linux R
load("H050_clone_info.RData")

clones = getTrees(clones, nproc=1,build="igphyml",exec="./software/igphyml/src/igphyml")
save(clones,file = "H050_clone_tree.RData")

load("./tree/H050_clone_tree.RData")

clones = getTrees(clones, nproc=1)

plots = plotTrees(clones, tips = "subType",tipsize = 6,palette=custom_palette)
plots[[1]]
treesToPDF(plots, file="./tree/Fig6D.H050_trees_subType.pdf", nrow=1, ncol=1,width = 4,height=3)

plots_time = plotTrees(clones, tips = "timepoint",tipsize = 6,palette=custom_palette)
plots_time[[1]]
treesToPDF(plots_time, file="./tree/Fig6D.H050_trees_timepoint.pdf", nrow=1, ncol=1,width = 4,height=3)

plots_ag = plotTrees(clones, tips = "antigen_main",tipsize = 6,palette=custom_palette)
plots_ag[[1]]
treesToPDF(plots_ag, file="./tree/Fig6D.H050_trees_antigen.pdf", nrow=1, ncol=1,width = 4,height=3)

plots_Iso = plotTrees(clones, tips = "Isotype",tipsize = 6,palette=custom_palette)
plots_Iso[[1]]
treesToPDF(plots_Iso, file="./tree/Fig6D.H050_trees_Isotype.pdf", nrow=1, ncol=1,width = 4,height=3)


# Prepare BCR data for Total HI
clone_info <- BC@meta.data[BC$clone_id %in% c(14124,6601) & BC$sample == "H052",]
clone_info$subType <- as.character(clone_info$subType)
clone_info$antigen_main <- as.character(clone_info$antigen_main)
names(clone_info) <- gsub("BCR_Heavy_", "", names(clone_info))

clone <- unique(clone_info$clone_id)

for (i in 1:length(clone)){
  print(clone[i])
  sub_clone_info <- clone_info[clone_info$clone_id == clone[i],]
  for (j in 1:nrow(sub_clone_info)){
    title <- sub_clone_info[j,"subType"]
    seq <- sub_clone_info[j,"germline_alignment"]
    if (is.character(seq) && length(seq) == 1) {
      seq <- strsplit(seq, "")[[1]]
    }
    write.fasta(sequences = list(seq), names = title, file.out = paste0("./tree/HI_germ/",clone[i],"_germ.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}


clone_info$germline_alignment_d_mask[clone_info$clone_id == 14124] <- "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCACTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAATATGTTTCAGCTATTAGTAGTAATGGGGGTAGCACATATTATGCAGACTCTGTGAAGGGCAGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTTCAAATGGGCAGCCTGAGAGCTGAGGACATGGCTGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clone_info$germline_alignment_d_mask[clone_info$clone_id == 6601] <- "CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clones = formatClones(clone_info, traits=c("Isotype","subType","timepoint","antigen_main"))
print(clones)

save(clones,file = "./tree/HI_clone_info.RData")

# use linux R
load("HI_clone_info.RData")
clones = getTrees(clones, nproc=1,build="igphyml",exec="./software/igphyml/src/igphyml")
save(clones,file = "HI_clone_tree.RData")

load("./tree/HI_clone_tree.RData")
clones = getTrees(clones, nproc=1)

plots = plotTrees(clones, tips = "subType",tipsize = 6,palette=custom_palette)
plots[[1]]
treesToPDF(plots, file="./tree/FigS7G.HI_trees_subType.pdf", nrow=1, ncol=2,width = 8,height=2.5)

plots_time = plotTrees(clones, tips = "timepoint",tipsize = 6,palette=custom_palette)
plots_time[[1]]
treesToPDF(plots_time, file="./tree/FigS7G.HI_trees_timepoint.pdf", nrow=1, ncol=2,width = 8,height=2.5)

plots_ag = plotTrees(clones, tips = "antigen_main",tipsize = 6,palette=custom_palette)
plots_ag[[1]]
treesToPDF(plots_ag, file="./tree/FigS7G.HI_trees_antigen.pdf", nrow=1, ncol=2,width = 8,height=2.5)

plots_Iso = plotTrees(clones, tips = "Isotype",tipsize = 6,palette=custom_palette)
plots_Iso[[1]]
treesToPDF(plots_Iso, file="./tree/FigS7G.HI_trees_Isotype.pdf", nrow=1, ncol=2,width = 8,height=2.5)


# Prepare BCR data for Total HVI
clone_info <- BC@meta.data[BC$clone_id %in% c(14666,14478,18650) & BC$sample %in% c("H002"),]
clone_info$subType <- as.character(clone_info$subType)
clone_info$antigen_main <- as.character(clone_info$antigen_main)
names(clone_info) <- gsub("BCR_Heavy_", "", names(clone_info))

clone <- unique(clone_info$clone_id)
for (i in 1:length(clone)){
  print(clone[i])
  sub_clone_info <- clone_info[clone_info$clone_id == clone[i],]
  for (j in 1:nrow(sub_clone_info)){
    title <- sub_clone_info[j,"subType"]
    seq <- sub_clone_info[j,"germline_alignment"]
    if (is.character(seq) && length(seq) == 1) {
      seq <- strsplit(seq, "")[[1]]
    }
    write.fasta(sequences = list(seq), names = title, file.out = paste0("./tree/HVI_germ/",clone[i],"_germ.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}
clone_info$germline_alignment_d_mask[clone_info$clone_id == 14478] <- "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGTAGCTATTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGGCCAACATAAAGCAAGATGGAAGTGAGAAATACTATGTGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clone_info$germline_alignment_d_mask[clone_info$clone_id == 14666] <- "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG"
clone_info$germline_alignment_d_mask[clone_info$clone_id == 18650] <- "GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"

clones = formatClones(clone_info, traits=c("Isotype","subType","timepoint","antigen_main"))
print(clones)

save(clones,file = "./tree/HVI_clone_info.RData")

# use linux R
load("HVI_clone_info.RData")
clones = getTrees(clones, nproc=1,build="igphyml",exec="./software/igphyml/src/igphyml")
save(clones,file = "HVI_clone_tree.RData")

load("./tree/HVI_clone_tree.RData")
clones = getTrees(clones, nproc=1)

plots = plotTrees(clones, tips = "subType",tipsize = 6,palette=custom_palette)
plots[[1]]
treesToPDF(plots, file="./tree/FigS7G.HVI_trees_subType.pdf", nrow=1, ncol=3,width = 12,height=2.5)

plots_time = plotTrees(clones, tips = "timepoint",tipsize = 6,palette=custom_palette)
plots_time[[1]]
treesToPDF(plots_time, file="./tree/FigS7G.HVI_trees_timepoint.pdf", nrow=1, ncol=3,width = 12,height=2.5)

plots_ag = plotTrees(clones, tips = "antigen_main",tipsize = 6,palette=custom_palette)
plots_ag[[1]]
treesToPDF(plots_ag, file="./tree/FigS7G.HVI_trees_antigen.pdf", nrow=1, ncol=3,width = 12,height=2.5)

plots_Iso = plotTrees(clones, tips = "Isotype",tipsize = 6,palette=custom_palette)
plots_Iso[[1]]
treesToPDF(plots_Iso, file="./tree/FigS7G.HVI_trees_Isotype.pdf", nrow=1, ncol=3,width = 12,height=2.5)


# Prepare BCR data for Total SLEI
clone_info <- BC@meta.data[BC$clone_id %in% c(8015) & BC$sample %in% c("R017"),]

clone_info$subType <- as.character(clone_info$subType)
clone_info$antigen_main <- as.character(clone_info$antigen_main)
names(clone_info) <- gsub("BCR_Heavy_", "", names(clone_info))

clone <- unique(clone_info$clone_id)

for (i in 1:length(clone)){
  print(clone[i])
  sub_clone_info <- clone_info[clone_info$clone_id == clone[i],]
  for (j in 1:nrow(sub_clone_info)){
    title <- sub_clone_info[j,"subType"]
    seq <- sub_clone_info[j,"germline_alignment"]
    if (is.character(seq) && length(seq) == 1) {
      seq <- strsplit(seq, "")[[1]]
    }
    write.fasta(sequences = list(seq), names = title, file.out = paste0("./tree/SLEI_germ/",clone[i],"_germ.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}

clone_info$germline_alignment_d_mask[clone_info$clone_id == 8015] <- "CAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTAGTAGTTACTACTGGGGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clones = formatClones(clone_info, traits=c("Isotype","subType","timepoint","antigen_main"))
print(clones)

save(clones,file = "./tree/SLEI_clone_info.RData")

# use linux R
load("SLEI_clone_info.RData")
clones = getTrees(clones, nproc=1,build="igphyml",exec="./software/igphyml/src/igphyml")
save(clones,file = "SLEI_clone_tree.RData")

load("./tree/SLEI_clone_tree.RData")
clones = getTrees(clones, nproc=1)

plots = plotTrees(clones, tips = "subType",tipsize = 6,palette=custom_palette)
plots[[1]]
treesToPDF(plots, file="./tree/FigS7G.SLEI_trees_subType.pdf", nrow=1, ncol=1,width = 4,height=2.5)

plots_time = plotTrees(clones, tips = "timepoint",tipsize = 6,palette=custom_palette)
plots_time[[1]]
treesToPDF(plots_time, file="./tree/FigS7G.SLEI_trees_timepoint.pdf", nrow=1, ncol=1,width = 4,height=2.5)

plots_ag = plotTrees(clones, tips = "antigen_main",tipsize = 6,palette=custom_palette)
plots_ag[[1]]
treesToPDF(plots_ag, file="./tree/FigS7G.SLEI_trees_antigen.pdf", nrow=1, ncol=1,width = 4,height=2.5)

plots_Iso = plotTrees(clones, tips = "Isotype",tipsize = 6,palette=custom_palette)
plots_Iso[[1]]
treesToPDF(plots_Iso, file="./tree/FigS7G.SLEI_trees_Isotype.pdf", nrow=1, ncol=1,width = 4,height=2.5)

# Prepare BCR data for Total SLEVI
clone_info <- BC@meta.data[BC$clone_id %in% c(4537,13353) & BC$sample %in% c("CR-3"),]
clone_info$subType <- as.character(clone_info$subType)
clone_info$antigen_main <- as.character(clone_info$antigen_main)
names(clone_info) <- gsub("BCR_Heavy_", "", names(clone_info))

clone <- unique(clone_info$clone_id)

for (i in 1:length(clone)){
  print(clone[i])
  sub_clone_info <- clone_info[clone_info$clone_id == clone[i],]
  for (j in 1:nrow(sub_clone_info)){
    title <- sub_clone_info[j,"subType"]
    seq <- sub_clone_info[j,"germline_alignment"]
    if (is.character(seq) && length(seq) == 1) {
      seq <- strsplit(seq, "")[[1]]
    }
    write.fasta(sequences = list(seq), names = title, file.out = paste0("./tree/SLEVI_germ/",clone[i],"_germ.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}

clone_info$germline_alignment_d_mask[clone_info$clone_id == 4537] <- "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clone_info$germline_alignment_d_mask[clone_info$clone_id == 13353] <- "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGTTATGAAATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGTGGTAGTACCATATACTACGCAGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
clones = formatClones(clone_info, traits=c("Isotype","subType","timepoint","antigen_main"))
print(clones)

save(clones,file = "./tree/SLEVI_clone_info.RData")

# use linux R
load("SLEVI_clone_info.RData")
clones = getTrees(clones, nproc=1,build="igphyml",exec="./software/igphyml/src/igphyml")
save(clones,file = "SLEVI_clone_tree.RData")

load("./tree/SLEVI_clone_tree.RData")
clones = getTrees(clones, nproc=1)

plots = plotTrees(clones, tips = "subType",tipsize = 6,palette=custom_palette)
plots[[3]]
treesToPDF(plots, file="./tree/FigS7G.SLEVI_trees_subType.pdf", nrow=1, ncol=3,width = 12,height=2.5)

plots_time = plotTrees(clones, tips = "timepoint",tipsize = 6,palette=custom_palette)
plots_time[[1]]
treesToPDF(plots_time, file="./tree/FigS7G.SLEVI_trees_timepoint.pdf", nrow=1, ncol=3,width = 12,height=2.5)

plots_ag = plotTrees(clones, tips = "antigen_main",tipsize = 6,palette=custom_palette)
plots_ag[[1]]
treesToPDF(plots_ag, file="./tree/FigS7G.SLEVI_trees_antigen.pdf", nrow=1, ncol=3,width = 12,height=2.5)

plots_Iso = plotTrees(clones, tips = "Isotype",tipsize = 6,palette=custom_palette)
plots_Iso[[1]]
treesToPDF(plots_Iso, file="./tree/FigS7G.SLEVI_trees_Isotype.pdf", nrow=1, ncol=3,width = 12,height=2.5)
