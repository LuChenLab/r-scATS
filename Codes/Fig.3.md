---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
suppressMessages({
  library("data.table")
  library("plyr")
  library("dplyr")
  library("ggplot2")
  library("Seurat")
  library(patchwork)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(scales)

})

```

```{r}
color3 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/color.celltype.Rds")

my_the <- theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_blank(),
        #axis.line = element_line(arrow = arrow(length = unit(0.5, 'cm'))),
        panel.background = element_rect(fill = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(fill=NA,color="black"),
        strip.text.x = element_text(size = 16),
        legend.position = "right")
pl <-  c("#FDAE61","#238B45")
names(pl) <- c("COVID19","Healthy")
color3 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/color.celltype.Rds")
cellorder <-  c("Class. mono","Inter. mono","Non-class. mono","pDC","mDC","NK","NKT","MAIT","CD8 memory","CD8 CTL","Treg","naïve CD4","naïve CD8","γδT", "Proliferative T","naïve B", "memory B", "plasma B" )


```


```{r}
scats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/00quant/DF.filtergene.Rds")
sub <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/00quant/NHC_NCov.Rds")

theta(scats) -> ats_psi
ats_psi[is.na(ats_psi)] <- 0
as(ats_psi, "sparseMatrix")->mtx

CreateAssayObject(mtx) ->sub[["ATS"]] #
DefaultAssay(sub) <- "ATS"
sub_ats <- ScaleData(sub, assay = 'ATS')

DefaultAssay(sub_ats) <- "RNA"

ElbowPlot(sub_ats) ->p
ggsave(p,file="/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN//elboe.pdf",height = 4,width = 7)
saveRDS(sub_ats,file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/DF.ATS.Rds")
saveRDS(sub_ats@assays$ATS@counts, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/raw.ATS.psi.Rds")
saveRDS(sub_ats@assays$RNA@data, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/RNA.data.Rds")

```


```{r,eval=F}
sub_ats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/DF.ATS.Rds")

np.rna = 10
np.ats <- 15
res.used <- 1

DefaultAssay(sub_ats) <- "ATS"
sub_ats <- RunPCA(sub_ats,features = rownames(sub_ats),reduction.name = 'ATSpsi.PCA',reduction.key = "ATSpsi PCA_")
sub_ats <- RunUMAP(sub_ats, dims = 1:np.ats, reduction = 'ATSpsi.PCA', assay = 'ATS', reduction.name = 'ATSpsi.UMAP', reduction.key = 'ATSpsi_UMAP_')

sub_ats <- FindMultiModalNeighbors(
  sub_ats, reduction.list = list("pca", "ATSpsi.PCA"), 
  dims.list = list(1:np.rna, 1:np.ats), 
  modality.weight.name = c("RNA.weight","ATSpsi.weight")#这里不同
)

sub_ats <- FindClusters(sub_ats, graph.name = "wsnn", algorithm = 3, resolution = res.used, verbose = FALSE)


sub_ats <- RunUMAP(sub_ats, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "WNN_UMAP_")
sub_ats <- RunUMAP(sub_ats, reduction = 'pca', dims = 1:np.rna, assay = 'RNA', reduction.name = 'RNA.UMAP', reduction.key = 'RNA_UMAP_')
sub_ats <- RunUMAP(sub_ats, reduction = 'ATSpsi.PCA', dims = 1:np.ats, assay = 'ATS', reduction.name = 'ATSpsi.UMAP2', reduction.key = 'ATSpsi_UMAP_')

saveRDS(sub_ats@meta.data,file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/wnn.meta.Rds")
saveRDS(sub_ats,file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/wnn.Rds")

```

# b-d

```{r}
my_the <- theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_blank(),
        #axis.line = element_line(arrow = arrow(length = unit(0.5, 'cm'))),
        panel.background = element_rect(fill = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(fill=NA,color="black"),
        strip.text.x = element_text(size = 16),
        legend.position = "none")

```


```{r}
color3 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/color.celltype.Rds")
library(ggpubr)
library(showtext)
DimPlot(sub_ats,label = T,repel = TRUE,label.size = 4,reduction = 'wnn.umap',group.by = "CellType")+ NoLegend()+
  scale_color_manual(values = c(color3))+
  my_the +
  guides(color = guide_legend( nrow = 2,byrow = TRUE,override.aes = list(size = 3)))+
  ggtitle("WNN") ->p1

DimPlot(sub_ats,label = T,repel = TRUE,label.size = 4,reduction = 'RNA.UMAP',group.by = "CellType")+ NoLegend()+
  scale_color_manual(values = c(color3))+
  my_the +
  guides(color = guide_legend( nrow = 2,byrow = TRUE,override.aes = list(size = 3)))+
  ggtitle("Gene expression") ->p2


DimPlot(sub_ats,label = T,repel = TRUE,label.size = 4,reduction = 'ATSpsi.UMAP2',group.by = "CellType")+ NoLegend()+
  scale_color_manual(values = c(color3))+
  my_the +
  guides(color = guide_legend( nrow = 2,byrow = TRUE,override.aes = list(size = 3)))+
  ggtitle("ATS PSI") ->p3

pdf( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/wnn.res3.pdf",height = 4,width = 12)
showtext_begin()
cowplot::plot_grid(p1,p2,p3,ncol = 3)

showtext_end()
dev.off()

```


# e
```{r, fig.width=4, fig.height=4}
meta <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/wnn.meta.Rds")
mean_weight <-
  data.frame(
    CellType = unique(meta$CellType),
    RNA = NA,
    ATS = NA
    
  )

mean_data <- lapply(mean_weight$CellType, function(x){
  subset_cell <- meta[meta$CellType %in% x,]
  data <- mean_weight[mean_weight$CellType %in% x,]
  data[,"RNA"] <- median(subset_cell$RNA.weight)
  data[,"ATS"] <- median(subset_cell$ATSpsi.weight)
  return(data)
})

mean_weight <- Reduce(rbind,mean_data)

mean_weight <- mean_weight[order(mean_weight$RNA, decreasing = T),]

rownames(mean_weight) <- mean_weight$CellType
mean_weight <- mean_weight[,-1]
mean_weight %>% tibble::rownames_to_column(.,var="cell") %>% melt(.) ->data_p

data_p
```


```{r, fig.width=3, fig.height=5}
library(showtext)
pdf( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/wnn.heatmap.pdf",height = 5,width = 3)
showtext_begin()
ggplot(data_p,aes(x=variable,y=factor(cell,levels = rownames(mean_weight)),fill=value))+
  geom_tile(color="grey", size=0.5, width=0.9, height=0.9)+
  scale_fill_gradient2(
    low = "#16499c", 
    high = "#e61f1b", 
    mid = "white",
    limits = c(0.3, 0.7),  # 手动设置范围（如c(0, 1)）
    midpoint = median(data_p$value)    # 或指定中点值
  )+
  theme_classic()+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank()
        )

showtext_end()
dev.off()

```

# f
```{r}
ats_marker <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/cell_DE.Rds")
ats_marker$p <- as.numeric(ats_marker$p)
ats_marker2 <- ats_marker[!is.na(ats_marker$p),]
ats_marker2$padj <- p.adjust(ats_marker2$p)
ats_marker3 <- ats_marker2[ats_marker2$padj< 0.05,]
ats_marker3[ats_marker3$group1=="naïve CD8",] ->df

df$ID <- mapIds(org.Hs.eg.db,
                  keys = df$gene,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals="first")

ora_res <- enrichGO(gene = df$ID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",#这里指定ID类型
                 ont = "BP", # "BP", "MF", "CC" 
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 5,# 最少的基因数量
                 maxGSSize = 300, # 最大的基因数量
                 readable = T # 把ENTREZID转换为SYMBOL
                 ) 

ora_res@result -> cd8_go
```

```{r fig.height=5, fig.width=8}
library(aPEAR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(cols4all)
set.seed(1)
p <- enrichmentNetwork(cd8_go[1:50,],
                        colorBy = 'p.adjust',
                        colorType = c("pval"),
                        nodeSize = 'Count',
                        fontSize = 4,
                        drawEllipses = T,
                        outerCutoff=0.3,
                       minClusterSize=4,
                        pCutoff = -50,
                        verbose = TRUE) +
  scale_color_continuous_c4a_div('spectral', reverse = F)

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/cell_DE_naive_cd8_go.pdf", height = 5, width = 8)
p
dev.off()

```


# h
```{r}
library(showtext)
source("/mnt/raid61/Personal_data/xuzijie/task/00scrna/01script/scrna_functions.R")

sub <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/wnn.Rds")

Idents(sub) <- sub$CellType
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/dotplot_SEMA4D.pdf", height = 3, width = 6.5)
showtext_begin()
i="SEMA4D"
my_rbind_plot(seu = sub, gene = i, ATS = c(rownames(sub@assays$ATS@counts)[grep(i, rownames(sub@assays$ATS@counts))]),ats_assay="ATS",
level = c(cellorder),lowcolor="lightgrey",highcolor="red")+coord_flip()

showtext_end()
dev.off()
```

# i
```{r}
sample_de_cell_all2 <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/sample_de_cell.Rds") %>% do.call(rbind,.)

sample_de_cell_all2[sample_de_cell_all2$type=="naïve CD8", ]-> samp_cd8

samp_cd8 %>% dplyr::group_by(condition) %>% dplyr::add_count(name="N")-> samp_cd8
samp_cd8$condition2 = paste0(samp_cd8$condition, " (", samp_cd8$N, ")")

label_data_left <- subset(samp_cd82, 
                          condition2 != "Common (15689)" & 
                          del_psi < 0 & 
                          gene %in% health_go2$gene)

label_data_right <- subset(samp_cd82, 
                           condition2 != "Common (15689)" & 
                           del_psi >= 0 & 
                           gene %in% covid_go2$gene)
mycol<- c( "#238B45", "#d9d9d9", "#FDAE61")
```


```{r fig.height=3.5, fig.width=4}
library(ggrepel)
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/sample_volcana2.pdf", height = 3.5, width = 4)

ggplot(samp_cd82, aes(x = del_psi, y = -log10(p), 
                color = factor(condition2, levels = c("Healthy (78)", "Common (15689)", "COVID19 (109)")))) +
  geom_point(alpha = 0.7,size=2) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7))+
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10)
  )+
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", linetype = "dashed", linewidth = 0.5) +
  ylim(0, 4)+geom_vline(xintercept = c(-0.1, 0.1), color = "darkgrey", linetype = "dashed", linewidth = 0.5) +
  geom_point(data = label_data_right,
             aes(x = del_psi, y = -log10(p)),
             color= '#e49d57', size = 4.5, alpha = 0.2) +
  geom_text_repel(data = label_data_right,
                  aes(x = del_psi, y = -log10(p), label = gene),
                  seed= 233,size= 3.5,
                  color= 'black',min.segment.length = 0,
                  force= 2,
                  force_pull= 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.5, #线段不透明度
                  nudge_x= 1 - label_data_right$del_psi, #标签x轴起始位置调整
                  direction= "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust= 0 #对齐标签：0右对齐，1左对齐，0.5居中
                  )+
  geom_point(data = label_data_left,
             aes(x = del_psi, y = -log10(p)),
             color= '#207d3e', size = 4.5, alpha = 0.2) +
  geom_text_repel(data = label_data_left,
                  aes(x = del_psi, y = -log10(p), label = gene),
                  seed= 233,size= 3.5,
                  color= 'black',min.segment.length = 0,
                  force= 2,
                  force_pull= 2,
                  box.padding = 0.1,
                  max.overlaps = 3,
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.5, #线段不透明度
                  nudge_x= -1 - label_data_left$del_psi, #标签x轴起始位置调整
                  direction= "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust= 1 #改为左对齐
                  )+
  labs(y="-Log10 (P-value)", x="(COVID-19 vs Healthy in naive CD8 cells)")
  
dev.off()

```

# j
```{r}
library(clusterProfiler)
library(org.Hs.eg.db)
samp_cd82 <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/DE_TSS.Rds")


mydf = data.frame("candigene" = samp_cd82[samp_cd82$condition!="Common",]$gene,
                  "type" = "DE") %>%
  .[!duplicated(.),]

mydf$ID <- mapIds(org.Hs.eg.db,
                  keys = mydf$candigene,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  
                  multiVals="first")

ora_res <- enrichGO(gene = mydf$ID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",#这里指定ID类型
                 ont = "BP", # "BP", "MF", "CC" 
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,# 最少的基因数量
                 maxGSSize = 300, # 最大的基因数量
                 readable = T # 把ENTREZID转换为SYMBOL
                 ) # 用什么函数富集都可以，主要是整理成 ora_res2 的格式

View(ora_res@result) 

saveRDS(ora_res, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/naive_cd8_go_all.Rds")


```
```{r}
go_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/naive_cd8_go_all.Rds")

data.frame(go_lst$type) -> plot2
plot2$labelx=rep(0,nrow(plot2))
plot2$labely=seq(nrow(plot2),1)
```


```{r fig.height=3, fig.width=4}

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/naive_cd8_go_all.pdf", height = 3, width = 4)
ggplot(data = plot2, 
       aes(-log10(pvalue), reorder(Description,-pvalue))) +
  geom_bar(stat="identity",
           # alpha=0.5,
           fill="#EED2EE",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = Description),
            size=3.5, 
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(qvalue)")+
  ggtitle("")+
  scale_x_continuous(expand = c(0,0))


dev.off()

```
