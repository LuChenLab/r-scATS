---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(scATS)
  library(GenomicRanges)
  library(magrittr)
  library(infercnv)
})

cellorder2 <- c("B","Mast","Epithilial")
color2 <- c("#aaddcd","#b6c2dd","#f2bddd")
names(color2) <- cellorder2
c("#7AC2E2", "#E56C83") -> co_cancer
order_cancer <- c("Non-malignant","Malignant")
names(co_cancer) <- order_cancer
```

# a

```{r}
seu <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/01malignant/new_seu.Rds")
gtfFile <- file.path("/mnt/raid64/ref_genomes/HomSap/release101/Homo_sapiens.GRCh38.101.sorted.gtf")
gtf <- rtracklayer::readGFFAsGRanges(gtfFile)

lapply(c("B"), function(cell){
  sub <- readRDS(paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/",cell, "/sub.Rds"))
  print(cell)
  mtx <- sub@assays$RNA@counts[rowSums(sub@assays$RNA@counts) >0 ,]
  
  gtf[gtf$type=="gene" & gtf$gene_name %in% rownames(sub@assays$RNA@counts),] %>% data.frame(.) %>% 
    select(seqnames,start,end,gene_name) %>% 
    .[!duplicated(.$gene_name),]->genes_df
  
  rownames(genes_df) <- genes_df$gene_name
  genes_df$gene_name <- NULL
  
  data.frame(sub@meta.data) %>% select(CellType_new_0.7) -> cellAnnotations
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = sub@assays$RNA@counts, # 基因原始矩阵
                                      annotations_file = cellAnnotations, # 细胞类型注释结果
                                      delim="\t",
                                      gene_order_file = genes_df,
                                      ref_group_names = cell)
  #注意：如果你没有参考细胞，可以设置 ref_group_names=NULL，在这种情况下，所有细胞基因的平均信号将用于定义基线（当做参照）
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/",cell,"/"),  # 如果输出文件夹不存在，会自己创建
                               cluster_by_groups=T,   # clusterx`
                               denoise=T,
                               HMM=T
                               )
  saveRDS(infercnv_obj,file = paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/", cell,"/", cell,".Rds"))
})

```

```{r}
seu <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/B_epi_seu.Rds")
cnv_lst_all <- read.table(paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/B/infercnv.observations.txt"))
colnames(cnv_lst_all) <- gsub("\\.","-", colnames(cnv_lst_all))
gtfFile <- file.path("/mnt/raid64/ref_genomes/HomSap/release101/Homo_sapiens.GRCh38.101.sorted.gtf")
gtf <- rtracklayer::readGFFAsGRanges(gtfFile)


# color
library(circlize)
col_fun <- colorRamp2(
  c(-4, 0, 4), 
  c("#323764", "white", "#69131A")
  )

col_fun(seq(-4, 4))


cellorder <- c("Non-malignant" ,"Malignant")
lapply(cellorder,function(cell){
  df0 <- data.frame("cellbc"=rownames(seu@meta.data)[seu@meta.data$infercnv2==cell],
                    "CellType"=cell)
  df0
})%>% do.call(rbind,.) -> df_col

df_col <- tibble::column_to_rownames(df_col,var = "cellbc") %>% data.frame(.)
df_col$CellType <- factor(df_col$CellType,levels = cellorder)

df_row <- data.frame(gtf[gtf$type == "gene",]) %>% dplyr::select(seqnames,gene_name) %>% .[!duplicated(.$gene_name),]
df_row[df_row$gene_name %in% rownames(cnv_lst_all),] -> df_row2
df_row2$gene_name -> rownames(df_row2)
df_row2$gene_name <- NULL

rowsplit <- structure(df_col$CellType,names=rownames(df_col))

# data
data2 <- cnv_lst_all[rownames(df_row2),rownames(df_col)]

# label
genes <- c(amp[amp %in% rownames(cnv_lst_all)],los[los %in% rownames(cnv_lst_all)])
label <- columnAnnotation(
  Zscore = anno_mark(at = which(rownames(data2) %in% genes),
                  labels = genes[genes %in% rownames(data2)],
                  labels_gp = gpar(fontsize = 10),
                  lines_gp = gpar())
  
)

```
```{r}
library(ComplexHeatmap)
df3 <- scale(t(data2),center = TRUE)
pdf( file ="/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/B_heatmap.all.pdf", height = 5, width = 9) 
showtext::showtext.begin()

Heatmap(df3,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        name = "infercnv",
        col = col_fun,
        na_col = "white",
        bottom_annotation = label,
        border = T,
        column_split = df_row2, # 只是分割开
        row_split = rowsplit)
showtext::showtext.end()
dev.off()

```

# b and c
```{r fig.width=3, fig.height=3}
CNV_score <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/B_cnvscore.Rds")
obj <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/B/17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj")
expr <- read.table(paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/","B","/infercnv.observations.txt"))

plot_chr_cluster(rownames(obj@gene_order[obj@gene_order$chr=="6",]), expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr6\ndeletion")


plot_chr_cluster(gene_lst[["11"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr11\namplification")-> chr11
plot_chr_cluster(gene_lst[["7"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr7\namplification")-> chr7
plot_chr_cluster(gene_lst[["20"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr20\namplification")-> chr20

plot_chr_cluster(gene_lst[["14"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr14\ndeletion")-> chr14
plot_chr_cluster(gene_lst[["6"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr6\ndeletion")-> chr6
plot_chr_cluster(gene_lst[["Y"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chrY\ndeletion")-> chrY
plot_chr_cluster(gene_lst[["21"]], expr, CNV_score[CNV_score$kmeans_class==9,"CB"], title = "chr21\ndeletion")-> chr21
```


```{r fig.width=4, fig.height=5}
library(ggpubr)
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V8/03infercnv/plot/B/chr_am_de.pdf", width = 4,height = 5 )
ggarrange(chr7,chr11, chr20,chr6, chr14, chr21, ncol = 3, nrow = 2, common.legend = T, legend = "top")

dev.off()
```
# d
```{r}
sub_ats <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/00quant/NSCLC.DF.ATS.Rds")
DefaultAssay(sub_ats) <- "RNA"
Idents(sub_ats) <- sub_ats@meta.data$infercnv2

markers <- FindAllMarkers(sub_ats, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(markers,file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/00quant/rna.marker.Rds")

```

```{r}
rna_marker <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/00quant/rna.marker.Rds")
rna_marker2 <- rna_marker[rna_marker$avg_log2FC >= 0.25 & rna_marker$p_val <= 0.05,]# 

ats_marker <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.DE_epi.Rds")

ats_marker$p <- as.numeric(ats_marker$p)
ats_marker$del.psi <- ats_marker$theta1 - ats_marker$theta2
ats_marker2 <- ats_marker[ats_marker$del.psi >= 0.1 & ats_marker$p<0.05,] #%>% dplyr::select(TSS, G1, Gene) %>% .[!duplicated(.),]
ats_marker2
```

### 2. overlap

```{r}
lapply(c("Non-malignant","Malignant"), function(i){
  both <- unique(intersect(rna_marker2[rna_marker2$cluster==i,]$gene, ats_marker2[ats_marker2$group1==i,]$gene))
  gene <- unique(setdiff(rna_marker2[rna_marker2$cluster==i,]$gene, ats_marker2[ats_marker2$group1==i,]$gene))
  ats <- unique(setdiff(ats_marker2[ats_marker2$group1==i,]$gene, rna_marker2[rna_marker2$cluster==i,]$gene))
  
  data.frame(gene = c(both, gene, ats),
             type = c(rep("ATS & RNA",length(both)), rep("RNA only",length(gene)), rep("ATS only",length(ats)))
             ) ->df
  data.frame(table(df$type)) %>% dplyr::mutate(cell = i,
                                               per = 100*round(.$Freq / sum(.$Freq),2))

}) %>% do.call(rbind,.)->markers
markers
```


```{r fig.height=4, fig.width=3.5}
library(dplyr)
markers$Var1 <- factor(markers$Var1, levels = c("RNA only", "ATS only", "ATS & RNA"))
cell_order <- markers %>%
  group_by(cell) %>%
  summarise(total = sum(Freq, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  pull(cell)

markers$cell <- factor(markers$cell, levels = cell_order)
```


```{r fig.height=4, fig.width=3.5}
pdf( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/rna.overlap.ats.pdf",height = 4,width = 3.5)
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
  "RNA only" = "#B7C9E2",   # 浅蓝灰
  "ATS only"  = "#2F8DBD",   # 蓝色
  "ATS & RNA" = "#0B6E4F"    # 深绿色
))+
  # scale_fill_brewer(palette = "Set3", label=c("Gene", "Overlap", "ATS"))+
  ylab("Genes numbers")+
  xlab(NULL)+
  theme(legend.position = "right",
        legend.title = element_blank(),
       axis.text.y = element_text(size = 13),
       axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
       legend.text = element_text(size = 13),
       axis.title  = element_text(size = 15))+
  guides(fill = guide_legend( nrow = 3))

showtext::showtext_end()
dev.off()
```

# e-g
```{r}
library(Seurat)
library(ggplot2)
library(ggpubr)

plot_ats_vln <- function(
  seu,
  feature,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "ATS",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "ATS value",
  show = c("p.format", "p.signif")[1]  # 显示 p 值 or 星号
) {
  df <- FetchData(seu, vars = c(feature, group.by), assay = assay)
  if(assay=="RNA"){feature=paste0("rna_",feature)}
  df <- df[!is.na(df[[feature]]), , drop = FALSE]
  if (length(unique(df[[group.by]])) < 2) {
    stop("Less than two groups after NA filtering. Cannot perform test.")
  }

  ymax <- max(df[[feature]], na.rm = TRUE)
  df$infercnv2 <- factor(df$infercnv2, levels = c("Non-malignant","Malignant"))
  p <- ggplot(df, aes(x = .data[[group.by]],
                      y = .data[[feature]],
                      fill = .data[[group.by]])) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    stat_compare_means(
      method = test.method,
      label = show,
      hide.ns = FALSE,
      label.y = ymax * 1.05,
      label.x = 1
    ) +
    scale_fill_manual(values = co_cancer)+
    theme_classic() +
    labs(x = NULL, y = ylab) +
    ggtitle(feature) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
    )

  return(p)
}
```

```{r fig.height=3, fig.width=4}
# sub_ats <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/00quant/NSCLC.DF.ATS.Rds")

ats="CCR6@6:167111796:+"
rna="CCR6"
plot_ats_vln(sub_ats,
  feature = ats,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "ATS",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "ATS value") -> p

plot_ats_vln(seu = sub_ats,
  feature = rna,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "RNA",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "RNA value") ->p2
ats_marker2[ats_marker2$TSS==ats,]
rna_marker2[rna_marker2$gene==rna,]

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/rna_ats_ccr6.pdf", height =3, width = 4.5)
cowplot::plot_grid(p2,p,ncol = 2)
dev.off()


ats="CCR2@3:46353861:+"
rna="CCR2"
plot_ats_vln(sub_ats,
  feature = ats,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "ATS",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "ATS value") -> p

plot_ats_vln(seu = sub_ats,
  feature = rna,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "RNA",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "RNA value") ->p2
ats_marker2[ats_marker2$TSS==ats,]
rna_marker2[rna_marker2$gene==rna,]

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/rna_ats_ccr2.pdf", height =3, width = 4.5)
cowplot::plot_grid(p2,p,ncol = 2)
dev.off()

```


```{r fig.height=3, fig.width=4}
# sub_ats <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/00quant/NSCLC.DF.ATS.Rds")

ats="CXCL16@17:4739902:-"
rna="CXCL16"
plot_ats_vln(sub_ats,
  feature = ats,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "ATS",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "ATS value") -> p

plot_ats_vln(seu = sub_ats,
  feature = rna,                # 基因名 or TSS 位点名（字符串）
  group.by = "infercnv2", # 分组列
  assay = "RNA",          # 默认 ATS assay
  test.method = "wilcox.test",  # 统计方法
  ylab = "RNA value") ->p2
ats_marker2[ats_marker2$TSS==ats,]
rna_marker2[rna_marker2$gene==rna,]
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/rna_only.pdf", height =3, width = 4.5)
cowplot::plot_grid(p2,p,ncol = 2)
dev.off()
```


# h
```{r}
###CSI
CSI.500bp.tss.mtx <- readRDS("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/luc/CSI.500bp.tss.mtx.Rds")
idx.2 <- apply(CSI.500bp.tss.mtx, 1, function(x){sum(x > 0.5) >= 0.05* ncol(CSI.500bp.tss.mtx)})

pdf("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/luc/CSI_all_TF_heatmap1.pdf", width = 8, height = 8)
p_csi <- pheatmap(CSI.500bp.tss.mtx[idx.2, idx.2],
                  fontsize_row = 16,
                  fontsize_col = 16,
                  show_rownames = T, 
                  show_colnames = T,
                  silent = T, 
                  cutree_rows = 4, 
                  cluster_cols = T,
                  cluster_rows = T,
                  color = colorRampPalette(c("#e5e4f3","white", "#9e0d0b"))(100))
print(p_csi)
dev.off()
```


# i and j
```{r}
tfs <- read.delim2("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/06TF/TFBSResult.tsv")

seu <- readRDS("/mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/NSCLC/nsclc_new_cluster.Rds") # new cluster 3
seu@assays$RNA@data -> mtx
count <- data.frame(sum=Matrix::rowSums(mtx))
```


```{r fig.height=4, fig.width=3}
data <- readRDS("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/luc/ELF1_gene.Rds")
data$count <- plyr::mapvalues(from = rownames(count), to = count$sum, x = data$target) %>% as.numeric()
ats <- c("RTKN2","CCR2","CCR6")

data$label <- ifelse(data$target %in% ats , 1,0)
data2 <- data[data$score < as.numeric(quantile(data$score, probs = 0.99)),]
data2$rank <- vctrs::vec_rank(-data2$score, ties="min")
data2$ID <- paste0("Rank:",data2$rank,"@",data2$target)


```


```{r fig.height=3, fig.width=3}
ggplot(data2, aes(x=score, y=log10(count)))+
  geom_point( linewidth = 1,size = 1) +
geom_point(data = data2[data2$label == "1",], 
           aes(x = score, y = log10(count), color="red"),size = 2)+
  theme_bw()+
  theme(axis.text = element_text(size=10),
      axis.title = element_text(size=12),
      legend.position = "none")+
  labs(x="Importance score of target genes",y="Gene expression (log10)",title = "ELF1")+
   geom_text_repel(data = data2[data2$label ==1,], 
                   aes(label = ID,color="red"),
                  size=3, max.overlaps = Inf, direction ="y",
                  hjust=45, vjust=3) +xlim(0,3.5)->p

```


```{r fig.height=4, fig.width=3}
data <- readRDS("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/luc/STAT1_gene.Rds")
data$count <- plyr::mapvalues(from = rownames(count), to = count$sum, x = data$target) %>% as.numeric()
ats <- c("RTKN2","CCR2","CCR6")

data$label <- ifelse(data$target %in% ats , 1,0)
data2 <- data[data$score < as.numeric(quantile(data$score, probs = 0.99)),]
data2$rank <- vctrs::vec_rank(-data2$score, ties="min")
data2$ID <- paste0("Rank:",data2$rank,"@",data2$target)

```
```{r fig.height=3, fig.width=3}
ggplot(data2, aes(x=score, y=log10(count)))+
  geom_point( linewidth = 1,size = 1) +
geom_point(data = data2[data2$label == "1",], 
           aes(x = score, y = log10(count), color="red"),size = 2)+
  theme_bw()+
  theme(axis.text = element_text(size=10),
      axis.title = element_text(size=12),
      legend.position = "none")+
  labs(x="Importance score of target genes",y="Gene expression (log10)",title = "STAT1")+
   geom_text_repel(data = data2[data2$label ==1,], 
                   aes(label = ID,color="red"),
                  size=3, max.overlaps = Inf, direction ="y",
                  hjust=45, vjust=3) +xlim(0,4.2)->p2


```


```{r fig.height=3, fig.width=6}
pdf("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/luc/elf1_stat1.pdf", height = 3,width = 6)
cowplot::plot_grid(p,p2, ncol = 2)
dev.off()
```

```{r}
data.frame(chr = c(rep("10",2), rep("3",3) ,rep("6",3)),
           start = c((62236121 -100 ), (62268845 -100), 
                     (46353861-500), (46354110-500), (46354125-500), 
                     (167111796-500), (167122753-500), (167123095-500)),
           end = c((62236121 +500 ), (62268845 +500), 
                     (46353861+100), (46354110+100), (46354125+100), 
                     (167111796+100), (167122753+100), (167123095+100)),
           TSS=c(candi)) -> df0

data.frame(chr = c(rep("10",2), rep("3",3) ,rep("6",3)),
           start = c((62236121 -1000 ), (62268845 -1000), 
                     (46353861-1000), (46354110-1000), (46354125-1000), 
                     (167111796-1000), (167122753-1000), (167123095-1000)),
           end = c((62236121 +1000 ), (62268845 +1000), 
                     (46353861+1000), (46354110+1000), (46354125+1000), 
                     (167111796+1000), (167122753+1000), (167123095+1000)),
           TSS=c(candi)) -> df

write.table(df, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/06TF/candi3_1000.bed", quote = F, sep = "\t", col.names = F, row.names = F)
```

```{r}
cd /mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/06TF/

bedtools getfasta -fi /mnt/raid66/Personal_data/xuzijie/ref/HomSap/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed candi3_1000.bed -fo candi3_1000.fasta

for i in MA0137.4 MA0473.4
do
echo $i
fimo --max-stored-scores 1000000 --oc ./ATS_out2/$i --thresh 1e-2 --alpha 1 --max-strand /mnt/raid61/Personal_data/xuzijie/ref/HomSap/JASPAR/MEME/$i.meme candi3_1000.fasta
done
```

```{r}
df0$name <- paste0(df0$chr, ":", df0$start, "-", df0$end)
df$name <- paste0(df$chr, ":", df$start, "-", df$end)
list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/06TF/ATS_out/", pattern = "fimo.txt", recursive = T, full.names = T) ->files
lapply(files, function(i){
  read.delim2(i)
}) %>% do.call(rbind,.) -> fimo_df
fimo_df$X.pattern.name <- ifelse(fimo_df$X.pattern.name=="MA0137.4","STAT1", "ELF1")
fimo_df$TSS <- plyr::mapvalues(from = df0$name, to=df0$TSS, x = fimo_df$sequence.name)
fimo_df$gene <- mapply(function(x) x[1],strsplit( fimo_df$TSS,"@"))
fimo_df[fimo_df$p.value <0.01,] ->fimo_df2
fimo_df2[fimo_df2$TSS=="CCR6@6:167123095:+" &fimo_df2$X.pattern.name=="STAT1",] %>% dplyr::select(start,matched.sequence,TSS)
```
```{r fig.width=5, fig.height=1.5}
library(ggseqlogo)

lapply(unique(fimo_df2$gene), function(gene){
  pdf(paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/06TF/",gene,".pdf"), height = 1.5, width = 3)
  ggseqlogo(fimo_df2[fimo_df2$gene==gene & fimo_df2$X.pattern.name=="STAT1",]$matched.sequence)->p1
  ggseqlogo(fimo_df2[fimo_df2$gene==gene & fimo_df2$X.pattern.name=="ELF1",]$matched.sequence)->p2
  cowplot::plot_grid(p1,p2,ncol = 2)
  dev.off()
})
```

