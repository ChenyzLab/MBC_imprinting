#### 00. Load packages and functions, define colors ####
library(openxlsx)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(gridExtra)

load("My_func.RData")

HI = "#839FBF"
HVI = "#074080"
SLEI = "#DFBFDF"
SLEVI = "#800080"
BF7 = "#88C0C0"
Cross = "grey90"
Ancestral = "grey95"
cluster_col <- c("A" = "#086233","B" = "#FF4040","C" = "#FFDF33","D1" = "#4EDDA2","D2" = "#4BAA30",
                 "E2" = "#AA364E","E3" = "#5C1A14","F1" = "#C58F50","F2" = "#325099","F3" = "#8AA3C8")

#### 01. UAMP :cluster####
mAb_embed <- read.csv("./DMS/processed_files/mAb_embed.csv")
sub_mAb_embed <- mAb_embed[mAb_embed$CohortType != "",]


highlight_antibodies <-  c('BRII-196','BD56-1854','REGN10987','BD57-1303',
                           'BD57-1271','BRII-198','BD55-5514', 'AZD1061')

p <- ggplot(mAb_embed, aes(x = UMAP1, y = UMAP2, color = manual_assign)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c(cluster_col), name = "Cluster") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.ticks = element_line(color = "black", size = 1),
    axis.text = element_text(size = 12)
  ) + 
  labs(title = "", x = "UMAP1", y = "UMAP2") + 
  coord_fixed(ratio = 0.95) + 
  geom_text_repel(
    data = subset(mAb_embed, antibody %in% highlight_antibodies),
    aes(label = antibody),
    fontface = "bold",
    size = 3,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "black",
    segment.size = 1,        
    arrow = arrow(length = unit(0.015, "npc"), type = "closed", ends = "last"),
    min.segment.length = 0,  
    force = 10,
    color = "black",
    max.overlaps = 20
  )

print(p)
ggsave("./DMS/figures/Fig.UMAP_mAb.pdf", plot = p, width = 5.5, height = 5)


#### 02. Visualize mAb info ####
#### 02.1 NT50 and RBD reactivity on UMAP ####
p <- ggplot(mAb_embed, aes(x = UMAP1, y = UMAP2,shape = BindingType,fill = BF7_NT50)) +
  geom_point(size = 4,stroke = 0) +
  scale_fill_gradient2(midpoint = 9000,high="grey", low="#F05060",name = "NT50 (ng/ml)") + 
  scale_shape_manual(
    values = c(25,24,23), 
    name = "RBD reactivity"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) + 
  labs(title = "", x = "", y = "") + 
  coord_fixed(ratio = 0.95)

print(p)
ggsave("./DMS/figures/Fig.UMAP_mAb_BF7NT50.pdf", plot = p, width = 5.5, height = 5)


p <- ggplot(mAb_embed, aes(x = UMAP1, y = UMAP2,shape = BindingType, fill = Ancestral_NT50)) +
  geom_point(size = 4,stroke = 0) +
  scale_fill_gradient2(midpoint = 9000,low="#F05060", high="grey",name = "NT50 (ng/ml)") + 
  scale_shape_manual(
    values = c(25,24,23), 
    name = "RBD reactivity"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) + 
  labs(title = "", x = "", y = "") +  
  coord_fixed(ratio = 0.95)

print(p)
ggsave("./DMS/figures/Fig.UMAP_mAb_AncestralNT50.pdf", plot = p, width =5.5, height = 5)


#### 02.2 Binder and Neutralizer distribution in each cluster ####
tt = as.data.frame(round(100*prop.table(x = table(mAb_embed$BindingType,mAb_embed$manual_main),margin = 2),2))
names(tt) = c("sample","cluster","Freq")
tt$sample = factor(tt$sample, levels = c("BF.7","Cross","Ancestral"))
p <- ggplot(tt,aes(y = Freq, x = factor(cluster), fill = sample)) + geom_bar(stat="identity")+
  scale_fill_manual(values = c(BF7,"grey87","grey92"),name = "RBD reactivity") +
  labs(x="", y="Percentage of RBD reactivity (%)") +
  theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,colour="black", size = rel(1.1)),axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(colour="black", size = rel(1.1)),legend.position = "right",axis.line.y = element_line(colour = "black", size = 0.5),
        legend.title = element_text(size = 10),panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA))
print(p)
ggsave("./DMS/figures/Fig.Barplot_BindingType_bycluster.pdf", plot = p, width = 4, height = 3)


tt = as.data.frame(round(100*prop.table(x = table(mAb_embed$Neutralizer,mAb_embed$manual_main),margin = 2),2))
names(tt) = c("sample","cluster","Freq")
tt$sample = factor(tt$sample, levels = c("BF.7","Cross","Ancestral","Non-Neutralizer"))

p <- ggplot(tt,aes(y = Freq, x = factor(cluster), fill = sample)) + geom_bar(stat="identity")+
  scale_fill_manual(values = c(BF7,"grey87","grey92","grey98"),name = "Neutralizer") +
  labs(x="", y="Percentage of neutralization (%)") +
  theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,colour="black", size = rel(1.1)),axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(colour="black", size = rel(1.1)),legend.position = "right",axis.line.y = element_line(colour = "black", size = 0.5),
        legend.title = element_text(size = 10),panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA))
print(p)
ggsave("./DMS/figures/Fig.Barplot_Neut_bycluster.pdf", plot = p, width = 4.3, height = 3)


#### 02.3 sample distribution and cluster distribution ####
sub_mAb_embed <- mAb_embed[mAb_embed$CohortType != "",]
sub_mAb_embed$CohortType <- factor(sub_mAb_embed$CohortType, levels = c("HI","HVI","SLEI","SLEVI"))

tt = as.data.frame(round(100*prop.table(x = table(sub_mAb_embed$manual_assign,sub_mAb_embed$CohortType),margin = 2),2))
names(tt) = c("cluster","sample","Freq")
tt$cluster <- factor(tt$cluster, levels = c("F2","F3","A","B","C","D2","D1","E2","E3","F1"))
p <- ggplot(tt,aes(y = Freq, x = sample, fill = cluster)) + geom_bar(stat="identity")+
  scale_fill_manual(values = c(cluster_col[c(9,10,1:3,5,4,6:8)]),name = "Cluster") +
  labs(x="", y="Percentage of epitope groups (%)") +
  theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,colour="black", size = rel(1.1)),axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(colour="black", size = rel(1.1)),legend.position = "right",axis.line.y = element_line(colour = "black", size = 0.5),
        legend.title = element_text(size = 10),panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA))
print(p)
ggsave("./DMS/figures/Fig.Barplot_cluster_bysample.pdf", plot = p, width = 4, height = 3)

#### 02.4 epitope composition of Cross and BF7 mAbs ####
sub_mAb_embed <- sub_mAb_embed[sub_mAb_embed$BindingType %in% c("Cross","BF.7"),]
sub_mAb_embed$sample_BindingType <- paste0(sub_mAb_embed$CohortType,"_",sub_mAb_embed$BindingType)
table(sub_mAb_embed$sample_BindingType)

tt = as.data.frame(round(100*prop.table(x = table(sub_mAb_embed$manual_assign,sub_mAb_embed$sample_BindingType),margin = 2),2))
names(tt) = c("cluster","sample","Freq")
tt$sample <- factor(tt$sample, levels = c("HI_Cross","HVI_Cross","SLEVI_Cross",
                                          "HI_BF.7","HVI_BF.7","SLEVI_BF.7"))
tt$cluster <- factor(tt$cluster, levels = c("F2","F3","A","B","C","D2","D1","E2","E3","F1"))
p <- ggplot(tt,aes(y = Freq, x = sample, fill = cluster)) + geom_bar(stat="identity")+
  scale_fill_manual(values = c(cluster_col[c(9,10,1:3,5,4,6:8)]),name = "Cluster") +
  labs(x="", y="Percentage of epitope groups (%)") +
  theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1,colour="black", size = rel(1.1)),axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(colour="black", size = rel(1.1)),legend.position = "right",axis.line.y = element_line(colour = "black", size = 0.5),
        legend.title = element_text(size = 10),panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA))
print(p)
ggsave("./DMS/figures/Fig.Barplot_epitope_distribution_inCrossBF7.pdf", plot = p, width = 4.5, height = 4)


#### 03. Plot mutation trend ####
dms_avg_cluster <- read.csv("./DMS/files/antibody_dms_merge_bycluster.csv")
data <- dms_avg_cluster %>%
  group_by(antibody, site) %>%
  mutate(mut_escape_site = sum(mut_escape)) %>%
  distinct(site, antibody, .keep_all = T) %>%
  dplyr::select(site, antibody, mut_escape_site)

p <- list()

data$site <- as.numeric(data$site)

trend_theme <- theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5, linetype = "dotted"),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.25, "cm"),
    
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, 
                              margin = margin(b = 20)),
    plot.margin = margin(t = 15, r = 20, b = 15, l = 15)
  )

cluster_col <- c("A" = "#086233","B" = "#FF4040","C" = "#FFD700","D1" = "#7CCD7C","D2" = "#4BAA30",
                 "E2" = "#AA364E","E3" = "#5C1A14","F1" = "#C58F50","F2" = "#325099","F3" = "#8AA3C8")

for (cluster in c("A", "B", "C", "D1", "D2", "E2", "E3", "F1", "F2", "F3")) {
  print(paste0("Plotting cluster ", cluster))
  
  data_cluster <- data %>% filter(antibody == cluster)
  top_sites <- data_cluster %>% 
    arrange(desc(mut_escape_site)) %>% 
    head(n = 8)
  
  print(min(top_sites$mut_escape_site))
  cluster_color <- cluster_col[cluster]
  p[[cluster]] <- ggplot(data_cluster, aes(x = site, y = mut_escape_site)) + 
    geom_ribbon(aes(ymin = 0, ymax = mut_escape_site), 
                fill = cluster_color, alpha = 0.15) +
    geom_line(color = cluster_color, size = 1.2, alpha = 0.9) + 
    geom_point(color = cluster_color, fill = "white", 
               shape = 21, size = 1.2, stroke = 1.2) + 
    geom_point(data = top_sites, 
               color = cluster_color, fill = cluster_color,
               shape = 21, size = 1.5, stroke = 1.5) +
    trend_theme +
    scale_x_continuous(breaks=seq(331,531,5)) + 
    scale_y_continuous(limits = c(0, max(data_cluster$mut_escape_site)+0.5)) +
    labs(
      y = 'Site Escape Score',
      x = 'RBD Residues'
    ) +
    geom_label_repel(
      data = top_sites, 
      aes(label = site), 
      min.segment.length = 0, 
      direction = "both", 
      fill = alpha("white", 0.9),
      color = "black",
      size = 3.5,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.5,
      max.overlaps = Inf
    )
  ggsave(paste0("./DMS/figures/FigS.Mut_Trend/", cluster, ".pdf"), plot = p[[cluster]], width = 8, height = 3)
}

#### 04. Plot Heatmap ####
library(circlize)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)

dms_avg <- read.csv("./DMS/files/antibody_dms_merge_avg.csv")
data <- dms_avg %>%
  group_by(antibody, site) %>%
  mutate(mut_escape_site = sum(mut_escape)) %>%
  distinct(site, antibody, .keep_all = T) %>%
  dplyr::select(site, antibody, mut_escape_site) %>%
  pivot_wider(
    id_cols = antibody,
    names_from = site,
    values_from = mut_escape_site
  )
data <- data.frame(data)
rownames(data) <- data$antibody

mAb_embed <- read.csv("./DMS/files/mAb_embed.csv")
mAb_embed <- mAb_embed[mAb_embed$CohortType != "",]
mAb_choose <- mAb_embed$antibody

data <- merge(data[data$antibody %in% mAb_choose,], mAb_embed[c(1,4:23)], by = "antibody")

colnames(data)[2:170] <- unique(paste0(dms_avg$wildtype,dms_avg$site))

#write.csv(data,"./DMS/figures/Fig4G.Heatmap/mAb_site_score_info.csv",row.names = F)
#data <- read.csv("./DMS/figures/Fig4G.Heatmap/mAb_site_score_info.csv")
data$antibody <- sapply(strsplit(as.character(data$antibody), "_"), `[`, 2)
data$BindingType <- factor(data$BindingType, levels = c("BF.7","Cross","Ancestral"))
data$CohortType <- factor(data$CohortType, levels = c("HI","HVI","SLEI","SLEVI"))
data <- data %>% arrange(manual_assign,BindingType,CohortType)
data$cluster <- data$manual_assign

escape_matrix <- data[,c(1:170)]
rownames(escape_matrix) <- escape_matrix$antibody
escape_matrix$antibody <- NULL


NT50_data <- data %>%
  dplyr::select(antibody, BF7_NT50, Ancestral_NT50) %>%
  distinct() %>%
  column_to_rownames("antibody")
colnames(NT50_data) <- c("BF7 NT50","Ancestral NT50")
NT50_data <- NT50_data[rownames(escape_matrix), ]

NT50_colors <- colorRamp2(
  breaks = c(min(NT50_data, na.rm = TRUE), 
             max(NT50_data, na.rm = TRUE)),
  colors = c("#F08080","#213558")
)

Binding_data <- data %>%
  dplyr::select(antibody, BindingType) %>%
  distinct() %>%
  column_to_rownames("antibody")

Binding_data <- Binding_data[rownames(escape_matrix), ]
Binding_colors <- c("BF.7" = BF7, "Cross" = "grey", "Ancestral" = "grey90")

Cohort_data <- data %>%
  dplyr::select(antibody, CohortType) %>%
  distinct() %>%
  column_to_rownames("antibody")

Cohort_data <- Cohort_data[rownames(escape_matrix), ]
Cohort_colors <- c("HI" = HI, "HVI" = HVI, "SLEI" = SLEI, "SLEVI" = SLEVI)

Cluster_data <- data %>%
  dplyr::select(antibody, cluster) %>%
  distinct() %>%
  column_to_rownames("antibody")

Cluster_data <- Cluster_data[rownames(escape_matrix), ]
Cluster_data <- factor(Cluster_data, levels = c("A", "B", "C", "D1", "D2", "E2", "E3", "F1", "F2", "F3"))

WT_ACE2_interface_site <- c("417", "446", "449", "453", "455", "456", "475", "486",
                            "487", "489", "493", "496", "498", "500", "501", "502", "505")

BF7_ACE2_interface_site <- c("403","417","449","453","455","456","475","476","477",
                             "486","487","489","490","493","494","498","500","501","502","505")

BF7_mutation_site <- c("339", "346", "371", "373", "375", "376","403", "405", "408", "417",
                       "440", "452", "477", "478", "484", "486","487","489", "498","500", "501","502","505")

# ACE2 interface
WT_ACE2_Interface <- data.frame(matrix("0", nrow = 1, ncol = ncol(escape_matrix)))
colnames(WT_ACE2_Interface) <- colnames(escape_matrix)
rownames(WT_ACE2_Interface) <- "WT_ACE2_Interface"

BF7_ACE2_Interface <- data.frame(matrix("0", nrow = 1, ncol = ncol(escape_matrix)))
colnames(BF7_ACE2_Interface) <- colnames(escape_matrix)
rownames(BF7_ACE2_Interface) <- "BF7_ACE2_Interface"

BF7_Mutation <- data.frame(matrix("0", nrow = 1, ncol = ncol(escape_matrix)))
colnames(BF7_Mutation) <- colnames(escape_matrix)
rownames(BF7_Mutation) <- "BF7_Mutation"

for (site in WT_ACE2_interface_site) {  
  WT_ACE2_Interface[, grepl(site, colnames(WT_ACE2_Interface))] <- "1"
}
for (site in BF7_ACE2_interface_site) {  
  BF7_ACE2_Interface[, grepl(site, colnames(BF7_ACE2_Interface))] <- "1"
}
for (site in BF7_mutation_site) {  
  BF7_Mutation[, grepl(site, colnames(BF7_Mutation))] <- "1"
}


WT_ACE2_Interface <- data.frame(t(WT_ACE2_Interface))
BF7_ACE2_Interface <- data.frame(t(BF7_ACE2_Interface))
BF7_Mutation <- data.frame(t(BF7_Mutation))

# HeatmapAnnotation
ACE2_anno <- HeatmapAnnotation(
  `Ancestral ACE2 Interface` = WT_ACE2_Interface$WT_ACE2_Interface,
  `BF7 ACE2 Interface` = BF7_ACE2_Interface$BF7_ACE2_Interface,
  `BF7 Mutation` = BF7_Mutation$BF7_Mutation,
  col = list(
    `Ancestral ACE2 Interface` = c("0" = "grey90", "1" = "#D90B29"), 
    `BF7 ACE2 Interface` = c("0" = "grey90", "1" = "#D90B29") , 
    `BF7 Mutation` = c("0" = "grey90", "1" = "black") 
  ),
  border = TRUE,
  gap = unit(1, "points"),
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10)
)


row_split = rep("A", nrow(escape_matrix))
row_split[5:6] = "B"
row_split[7] = "C"
row_split[8:12] = "D1"
row_split[13:15] = "D2"
row_split[16:22] = "E2"
row_split[23:27] = "E3"
row_split[28:31] = "F1"
row_split[32:39] = "F2"
row_split[40] = "F3"
row_split <- factor(row_split, levels = c("A","B","C","D1","D2","E2","E3","F1","F2","F3"))


main_heatmap <- Heatmap(as.matrix(escape_matrix),
                        name = "Escape Score",
                        col = colorRamp2(seq(from=0,to=8,length=5),rev(brewer.pal(5, "Spectral"))),
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        row_title = NULL,
                        cluster_rows = F,
                        cluster_columns = F,
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10),
                        row_split = row_split,
                        border = T,
                        row_names_side ="left",
                        column_names_side = "top",
                        use_raster=F,
                        column_gap = unit(1.5, "mm"),
                        width = unit(65,"cm"),
                        height = unit(20,"cm"),
                        row_labels = rownames(escape_matrix),
                        row_title_rot = 0,
                        column_names_rot = 90,
                        bottom_annotation = ACE2_anno,
)

Cluster_heatmap <- Heatmap(
  as.matrix(Cluster_data),
  name = "Cluster",
  col = cluster_col,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  width = unit(0.7, "cm"),
  height = unit(20, "cm"),
  row_split = row_split,
  border = T,
  row_names_side ="left",
  use_raster=F,
  column_gap = unit(1.5, "mm"),
  row_labels = rownames(escape_matrix),
  row_title_rot = 0,
  show_heatmap_legend = FALSE
)

NT50_heatmap <- Heatmap(
  as.matrix(NT50_data),
  name = "NT50",
  col = NT50_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  width = unit(1.4, "cm"),
  height = unit(20, "cm"),
  column_names_side = "top"
)

Binding_heatmap <- Heatmap(
  as.matrix(Binding_data),
  name = "RBD reactivity",
  col = Binding_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(0.7, "cm"),
  height = unit(20, "cm")
)

Cohort_heatmap <- Heatmap(
  as.matrix(Cohort_data),
  name = "CohortType",
  col = Cohort_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(0.7, "cm"),
  height = unit(20, "cm")
)

combined_heatmap <- Cluster_heatmap + main_heatmap + NT50_heatmap + Binding_heatmap + Cohort_heatmap

pdf("./DMS/figures/Fig.Heatmap/Heatmap_mAb_site_score.pdf",width = 30,height = 12)
draw(combined_heatmap, heatmap_legend_side = "right")
dev.off()


