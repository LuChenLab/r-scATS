---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
suppressMessages({
  library(scATS)
  library(Seurat)
  library(ggplot2)
  library(parallel)
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(magrittr)  
  library(ggpubr)
})
```


```{r}
color3 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/color.celltype.Rds")
# write.table(data.frame(color3), file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/color_celltype.txt", quote = F, col.names = F, row.names = T, sep = "\t")
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

# a
```{r}
sub <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/00quant/NHC_NCov.Rds")
Idents(sub) <- sub@meta.data$CellType

```

```{r}

`naïve B`= c("IGHM","TCL1A")
`memory B`= c("MS4A1","CD79A")
`plasma B`= c("MZB1","JCHAIN")
`naïve CD4`=c("CD4","IL7R")
`naïve CD8`= c("CD8A","CCR7","CD8B")
Treg= c("SELL","IL32","FOXP3")
MAIT = c("TRAV1-2","GZMK")
γδT = c("TRGV9", "TRDV2")
`Proliferative T`= c("MKI67","STMN1")
`CD8 memory`= c( "GZMB","SYNE1")
`CD8 CTL` = c( "GZMA")
NKT = c("NKG7","CCL4")
NK = c("GNLY","SPON2")
`Class. mono` = c("LYZ","CD14")
`Inter. mono`= c("IFITM3")
`Non-class. mono`= c("FCGR3A")
mDC= c("CD1C")
pDC= c("CLEC4C")

  
marker <- list(
  "naïve B"= `naïve B`,
  "memory B"= `memory B`,
  "plasma B"= `plasma B`,
  "naïve CD4"= `naïve CD4`,
  "naïve CD8"= `naïve CD8`,
  "Treg"= Treg,
  "MAIT" = MAIT,
  "γδT" = γδT,
  "Proliferative T"=`Proliferative T`,
  "CD8 memory"= `CD8 memory`,
  "CD8 CTL" = `CD8 CTL`,
  "NKT" = NKT,
  "NK" = NK,
  "Class. mono" = `Class. mono`,
  "Inter. mono"= `Inter. mono`,
  "Non-class. mono"= `Non-class. mono`,
  "mDC"= mDC,
  "pDC"= pDC
)

DefaultAssay(sub) <- "RNA"
sub@meta.data$CellType <- factor(sub@meta.data$CellType,levels = levels(df$cluster))

```


```{r fig.height=10,fig.width=7}
source("/mnt/raid61/Personal_data/xuzijie/task/00scrna/01script/scrna_functions.R")
sub@meta.data$CellType <- factor(sub@meta.data$CellType,levels = names(marker))

Idents(sub) <- sub$CellType
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/00quant/dotplot.pdf", height = 10, width = 7)
showtext::showtext.begin()
my_dotplot(marker, sub)
showtext::showtext.end()
dev.off()
```


# b
```{r fig.height=3, fig.width=6}

Idents(sub) <- sub$group
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/group_umap.pdf", height = 3, width = 6)
DimPlot(sub,label = F,repel = TRUE,label.size = 4,reduction = 'wnn.umap',split.by = "group")+ NoLegend()+
  scale_color_manual(values = c(pl))+
  my_the +
  guides(color = guide_legend( nrow = 2,byrow = TRUE,override.aes = list(size = 3)))
dev.off()
```

# c
```{r}
rna_marker <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V3/01DE/cell_DE/rna.marker.Rds")
rna_marker2 <- rna_marker[rna_marker$avg_logFC >= 0.25 & rna_marker$p_val <= 0.05,]# 


ats_marker <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/cell_DE.Rds") 
ats_marker$p <- as.numeric(ats_marker$p)

ats_marker$del.psi <- ats_marker$theta1 - ats_marker$theta2
ats_marker2 <- ats_marker[ats_marker$del.psi >= 0.1 & ats_marker$p <= 0.05 ,] #%>% dplyr::select(TSS, G1, Gene) %>% .[!duplicated(.),]
ats_marker2
```

```{r}
lapply(unique(ats_marker2$group1), function(i){
  both <- unique(intersect(rna_marker2[rna_marker2$cluster==i,]$gene, ats_marker2[ats_marker2$group1==i,]$Gene))
  gene <- unique(setdiff(rna_marker2[rna_marker2$cluster==i,]$gene, ats_marker2[ats_marker2$group1==i,]$Gene))
  ats <- unique(setdiff(ats_marker2[ats_marker2$group1==i,]$Gene, rna_marker2[rna_marker2$cluster==i,]$gene))
  
  data.frame(gene = c(both, gene, ats),
             type = c(rep("Overlap",length(both)), rep("Gene_Expression",length(gene)), rep("ATS_PSI",length(ats)))
             ) ->df
  data.frame(table(df$type)) %>% dplyr::mutate(cell = i,
                                               per = 100*round(.$Freq / sum(.$Freq),2))

}) %>% do.call(rbind,.)->markers

markers[markers$Var1=="Gene_Expression",]->rna_n
rna_n[order(rna_n$per, decreasing = T),]$cell ->ord    


```


```{r fig.height=4, fig.width=6}
library(dplyr)
markers$Var1 <- factor(markers$Var1, levels = c("Gene_Expression", "ATS_PSI", "Overlap"))
cell_order <- markers %>%
  group_by(cell) %>%
  summarise(total = sum(Freq, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  pull(cell)

markers$cell <- factor(markers$cell, levels = cell_order)


pdf( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/01DE/cell_DE/rna.overlap.ats2.pdf",height = 4,width = 7)
showtext::showtext_begin()
ggplot(markers,
       aes(x=factor(cell, cell_order),
           y=Freq,
           fill=Var1
           ))+
geom_bar(stat = "identity",size=0.2, width = 0.9)+
  theme_classic()+
  #scale_fill_manual(values = c(color))+
  scale_fill_manual(values = c(
  "Gene_Expression" = "#B7C9E2",   # 浅蓝灰
  "ATS_PSI"  = "#2F8DBD",   # 蓝色
  "Overlap" = "#0B6E4F"    # 深绿色
))+
  # scale_fill_brewer(palette = "Set3", label=c("Gene", "Overlap", "ATS"))+
  ylab("Genes numbers")+
  xlab(NULL)+
 theme(legend.position = "top",
        legend.title = element_blank(),
       axis.text.y = element_text(size = 12),
       axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
       legend.text = element_text(size = 12),
       axis.title  = element_text(size = 14))+
  guides(fill = guide_legend( nrow = 1))

showtext::showtext_end()
dev.off()
```


# e and f

```{bash}
pyscenic grn --num_workers 15 \
--sparse \
--method grnboost2 \
--output sce.adj.csv \
sce.loom \
/mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/SCENIC/hs_hgnc_tfs.txt # tfs 


pyscenic ctx --num_workers 15 \
--output sce.regulons.csv \
--expression_mtx_fname sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
sce.adj.csv \
/mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


pyscenic aucell --num_workers 40 \
--output sce_SCENIC_xzj.loom \
sce.loom \
sce.regulons.csv

```

```{r}
library(SCENIC)
meta <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V3/00quant/DF_UTR.colData.Rds") %>% data.frame(.)

cellanno_lineage <- meta %>% [meta$CellType=="naïve CD8",] %>%
  dplyr::select(group)

sce_SCENIC <- open_loom("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/covid19/sce_SCENIC2.loom")

regulons.500bp.tss <- get_regulons(sce_SCENIC, column.attr.name="Regulons") # 行为regulon，列为基因的ont-hot矩阵，行的加和即为regulonsToGeneLists的结果

regulonsGene.500bp.tss <- regulonsToGeneLists(regulons.500bp.tss)#每个regulon靶向的基因
regulonAUC.500bp.tss <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC') # 每个cellbc的regulon活性


regulonAUC.500bp.tss2 <- regulonAUC.500bp.tss[,rownames(cellanno_lineage)]

RSS.500bp.tss.lineage <- SCENIC::calcRSS(AUC = AUCell::getAUC(regulonAUC.500bp.tss2),cellAnnotation = cellanno_lineage[colnames(AUCell::getAUC(regulonAUC.500bp.tss2)),"group"])

```

```{r fig.height=3, fig.width=2.5}
RSS.500bp.tss.lineage
rssPlot <- SCENIC::plotRSS(
  RSS.500bp.tss.lineage,
  zThreshold = 0,
  cluster_columns = FALSE,
  order_rows = TRUE,
  thr=0,
  varName = "group"
)

rssPlot$plot
```

```{r}
sub2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/04WNN/naive_cd8.Rds")
Idents(sub2) <- sub2$group
DotPlot(sub2, features = c("TCF7"), assay = "RNA")->p0
```

```{r fig.height=3, fig.width=7}
rssPlot$plot +coord_flip() ->p2
p0+theme_bw() -> p1
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/06TF/covid_tcf7.pdf", height = 3, width = 7)
cowplot::plot_grid(p1,p2, nrow = 1, rel_widths = c(0.5,0.5))
dev.off()
```

