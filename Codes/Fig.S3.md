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

# a

```{r}
MOFAobject <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/MOFA2_bm_theta_0424.Rds")
seu <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/bm.DF.ATS.Rds")
seu@meta.data %>% tibble::rownames_to_column(., var = "sample") -> bm_meta

merge(MOFAobject@samples_metadata, bm_meta) ->df

MOFAobject@samples_metadata <- df

MOFAobject@samples_metadata$Lineages <- ifelse(MOFAobject@samples_metadata$wsnn_cell_type %in% c("CMP", "BasophilP", "HSC", "GMP", "MKP", "EBP"), "Progenitor", "Mature")
MOFAobject@samples_metadata$Lineages <- factor(MOFAobject@samples_metadata$Lineages, c("Progenitor", "Mature"))
saveRDS(MOFAobject,"/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/MOFA2_bm_theta_0424.Rds")

```

```{r fig.height=4, fig.width=4}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/bm.cor.factors_theta.pdf",height = 4, width = 4)
plot_factor_cor(MOFAobject)
dev.off()
```

# b
```{r fig.height=3, fig.width=6}
sub <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03WNN/wnn.Rds")
Idents(sub) <- sub$wsnn_cell_type

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/dotplot_Capg.pdf", height = 3, width = 6)

gene <- "Capg"
my_rbind_plot(seu = sub, gene = gene, ATS = c(rownames(sub@assays$ats_Theta@counts)[grep(gene, rownames(sub@assays$ats_Theta@counts))])[1:2],
              ats_assay = "ats_Theta", level = c(cellorder),lowcolor="lightgrey",highcolor="red") +coord_flip()

dev.off()

```


# c
```{r}
gtfFile <- "/mnt/raid62/Tabula_Muris/mm10_fa/Mus_musculus.GRCm38.93.gtf"
txdb_smart <- makeTxDbFromGFF(gtfFile)
gr <- range(unique(exons(txdb_smart)[queryHits(findOverlaps(exons(txdb_smart), as("6:72544467-72555272:+", "GRanges"), maxgap = 1))]))
gr <- keepSeqlevels(gr, as.character(runValue(seqnames(gr))))

```

```{r}
gtf <- data.table::as.data.table(rtracklayer::readGFF(gtfFile))
gr.txdb_smart <- biovizBase::crunch(txdb_smart, which = gr)
gr.txdb_smart <- split(gr.txdb_smart, gr.txdb_smart$tx_name)
names(gr.txdb_smart) <- plyr::mapvalues(names(gr.txdb_smart), gtf[type == "transcript", transcript_id], gtf[type == "transcript", transcript_name], warn_missing = F)
sk <- mapply(function(x) {
  y <- subset(x, type == "exon")
  paste(as.character(sort(ranges(y))), collapse = ",")
}, gr.txdb_smart)
sk <- split(names(sk), sk)
names(sk) <- mapply(function(x) paste(sort(x), collapse = ", "), sk)

gr.txdb_smart <- lapply(sk, function(x) gr.txdb_smart[[which(names(gr.txdb_smart) == x[1])]])
gr.txdb_smart <- GRangesList(gr.txdb_smart)
gr.txdb_smart <- gr.txdb_smart[!names(gr.txdb_smart) %in% c("Ptprc-204", "Ptprc-205")]

p8 <- ggbio::autoplot(gr.txdb_smart)@ggplot + ggbio::theme_null()
p8
```

```{r}
CellMeta <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/candi/TC/cellmeta_tablus.Rds")
icas <- readRDS("/mnt/raid61/Personal_data/tangchao/40TCellsQTL5/analysis/12.SingleCell/03.Tabula_Muris_Marrow_HSPC/RData/ICAS_Object.Rds")

icas@ASType[which(icas@HostGene$gene_name=="Capg"),]
as.character(icas@HostGene[which(icas@HostGene$gene_name=="Capg"),]) -> capg

names(icas@ASType) <- as.character(icas@ASType)
icas@ASType[capg,]
CellMeta <- CellMeta[file.exists(paste0("/mnt/raid62/Tabula_Muris/STAR/SmartSeq2/Marrow/", Run, "Aligned.sortedByCoord.out.bam")), ]
cells <- CellMeta[, unique(Cell)]
```

```{r}
gr <- GRanges(seqnames(gr), IRanges(min(layer_scales(p8)$x$range$range), max(layer_scales(p8)$x$range$range)), strand = strand(gr))
library(parallel)
minMapQuality <- 0
param0 = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isNotPassingQualityControls = FALSE), 
                                 what = c("qname", "seq", "mapq", "qual"),
                                 which = gr, 
                                 mapqFilter = minMapQuality)
CovList <- mclapply(cells, function(i) {
  run <- CellMeta[Cell == i, Run]
  if(length(run) == 1) {
    fils <- paste0("/mnt/raid62/Tabula_Muris/STAR/SmartSeq2/Marrow/", run, "Aligned.sortedByCoord.out.bam")
    fils <- fils[file.exists(fils)]
    bam <- tryCatch(GenomicAlignments::readGAlignments(file = fils, param = param0), error = function(e) NULL)
    if(is.null(bam)) {
      rep(0, length(start(gr):end(gr)))
    } else {
      cov <- as.integer(coverage(bam)[[seqlevels(gr)]][start(gr):end(gr)])
    }
  } else {
    fils <- paste0("/mnt/raid62/Tabula_Muris/STAR/SmartSeq2/Marrow/", run, "Aligned.sortedByCoord.out.bam")
    fils <- fils[file.exists(fils)]
    cov <- lapply(fils, function(x) {
      bam <- tryCatch(GenomicAlignments::readGAlignments(file = x, param = param0), error = function(e) NULL)
      if(is.null(bam)) {
        rep(0, length(start(gr):end(gr)))
      } else {
        as.integer(coverage(bam)[[seqlevels(gr)]][start(gr):end(gr)])
      }
    })
    cov <- colSums(do.call(rbind, cov))
  }
  return(cov)
}, mc.cores = 40)
Cov <- as.data.frame(do.call(rbind, CovList))
row.names(Cov) <- cells
colnames(Cov) <- start(gr):end(gr)


```

```{r}
melt_Cov <- data.table::melt.data.table(data.table::as.data.table(Cov, keep.rownames = "Cell"), variable.name = "Pos", value.name = "Depth")
melt_Cov[, Pos := as.numeric(as.character(Pos))]
melt_Cov <- melt_Cov[, .(Pos = Pos, Depth = Depth, Coverage = Depth/max(Depth)), by = Cell]
melt_Cov[is.na(Coverage), Coverage := 0]

melt_Cov <- merge(melt_Cov, unique(CellMeta[, .(Cell, cell)]), by = "Cell")

melt_Cov[, cell := factor(cell, levels = levels(SummarizedExperiment::colData(icas)[[icas@design]]))]
data.table::setkey(melt_Cov, cell)
Cell_od <- unique(melt_Cov[, .(Cell, cell)])[, Cell]
melt_Cov[, Cell := factor(Cell, levels = Cell_od)]
melt_Cov[, cell := factor(cell, levels = levels(SummarizedExperiment::colData(icas)[[icas@design]]))]

C2U <- melt_Cov[, max(Coverage), Cell][V1 > 0, as.character(Cell)]
```


```{r}
anncolor <- list(Cell = RColorBrewer::brewer.pal(n = length(levels(SummarizedExperiment::colData(icas)[[icas@design]])), name = "Paired"))
names(anncolor$Cell) <- levels(SummarizedExperiment::colData(icas)[[icas@design]])
ggplot(data = melt_Cov, 
       mapping = aes(x = Pos, y = Cell, fill = Coverage)) + 
  geom_tile(show.legend = FALSE) + 
  # theme_bw(base_size = 15) +
  theme_classic(base_size = 15) +
  scale_fill_gradient(high = "#EF3B2C", low = "#FFF5F0", na.value = "white") + 
  ggforce::facet_col(vars(cell), scales = 'free_y', space = 'free', strip.position = "left") + 
  scale_x_continuous(expand = c(0.001, 0.001)) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line.x = element_blank(), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 12), 
        strip.background = element_rect(fill = anncolor$Cell, colour = NA), 
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0)) -> p6

g <- ggplot_gtable(ggplot_build(p6))
stripl <- which(grepl('strip-l', g$layout$name))

fills <- anncolor$Cell
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k + 1
}
p6 <- as_ggplot(g)

```

```{r}
yls <- layer_scales(p8)$y$range$range
p8 <- p8 + scale_y_continuous(limits = c(min(yls), max(yls) + 0.5))
p8 <- p8 + theme(plot.margin = unit(c(0, 0.5, 1, 6.1), "line"))
pe <- cowplot::plot_grid(p6, p8, ncol = 1, rel_heights = c(6, 3.5))

```


```{r}
ggplot(data = melt_Cov2, 
       mapping = aes(x = Pos, y = Cell, fill = Coverage)) + 
  geom_tile(show.legend = FALSE) + 
  # theme_bw(base_size = 15) +
  theme_classic(base_size = 15) +
  scale_fill_gradient(high = "#EF3B2C", low = "#FFF5F0", na.value = "white") + 
  ggforce::facet_col(vars(cell), scales = 'free_y', space = 'free', strip.position = "left") + 
  scale_x_continuous(expand = c(0.001, 0.001)) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line.x = element_blank(), 
        panel.grid = element_blank(), 
        strip.text = element_text(size = 12), 
        strip.background = element_rect(fill = anncolor$Cell, colour = NA), 
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0)) -> p6_2

g_2 <- ggplot_gtable(ggplot_build(p6_2))
stripl <- which(grepl('strip-l', g_2$layout$name))

fills <- anncolor$Cell
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g_2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k + 1
}
p6_2 <- as_ggplot(g_2)
p6_2

```

```{r}
Mat2 <- melt_Cov[Cell %in% C2U, .(me = mean(Coverage, na.rm = TRUE), v = sd(Coverage, na.rm = TRUE)), by = list(cell, Pos)]
setkey(Mat2, cell, Pos)
Mat2 <- merge(Mat2, do.call(rbind, lapply(Mat2[, levels(cell)], function(x) data.table::data.table(cell = x, Pos = melt_Cov[, min(Pos)]:melt_Cov[, max(Pos)]))), by = c("cell", "Pos"), all.y = TRUE)
Mat2[is.na(me), me := 0]
Mat2[is.na(v), v := 0]

ggplot(Mat2, aes(x = Pos, y = me)) +
  geom_step(aes(colour = cell), fill = NA) +
  # geom_ribbon(aes(x = Pos, ymin = me - v, ymax = me + v, fill = cell), alpha = 0.2) +
  scale_x_continuous(expand = c(0.001, 0.001)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        legend.position = "none") +
  labs(y = "Scaled\nexpression") + 
  scale_color_manual(values = anncolor$Cell) -> p7
p7
```

```{r fig.width=14, fig.height=1.5}
yls <- layer_scales(p8)$y$range$range
p8 <- p8 + scale_y_continuous(limits = c(min(yls), max(yls) + 0.5))
```

```{r fig.width=14, fig.height=6}
p6 <- p6 + theme(plot.margin = unit(c(1, 0, 0, 1), "line"))
p7 <- p7 + theme(plot.margin = unit(c(0, 0.5, 0, 4.21), "line"))
cowplot::plot_grid(p6, p7, ncol = 1, rel_heights = c(6, 1.5), labels = c("e"), label_size = 30, label_y = 1.04)
```

```{r fig.width=14, fig.height=8}
p8 <- p8 + scale_x_continuous(expand = c(0, 0))
p8 <- p8 + theme(plot.margin = unit(c(0, 0.5, 1, 8.5), "line"))
pe <- cowplot::plot_grid(p6, p7, p8, ncol = 1, rel_heights = c(6, 1.5, 1.75), label_size = 30, label_y = 1.04)
pe
```

```{r}
bams <- mclapply(levels(SummarizedExperiment::colData(icas)[[icas@design]]), function(i) {
  run <- CellMeta[cell == i, Run]
  fils <- paste0("/mnt/raid62/Tabula_Muris/STAR/SmartSeq2/Marrow/", run, "Aligned.sortedByCoord.out.bam")
  bams <- lapply(fils, function(x) {
    tryCatch(GenomicAlignments::readGAlignments(file = x, param = param0), error = function(e) NULL)
  })
  do.call(c, bams)
}, mc.cores = 40)
names(bams) <- levels(SummarizedExperiment::colData(icas)[[icas@design]])
```

```{r}
CovsList <- mclapply(bams, function(bam) {
  data.table::data.table(x = start(gr):end(gr), y = as.integer(coverage(bam)[[seqlevels(gr)]][start(gr):end(gr)]))
}, mc.cores = 40)
```

```{r fig.width=14, fig.height=4}
CovsTab <- data.table::data.table(cell = rep(names(CovsList), mapply(nrow, CovsList)), do.call(rbind, CovsList))
CovsTab[, cell := factor(cell, levels = names(anncolor$Cell))]
ggplot(CovsTab, aes(x, y, fill = cell)) + 
  geom_area() + 
  scale_fill_manual(values = anncolor$Cell, guide = "none") + 
  ggforce::facet_col(vars(cell), scales = 'free_y', strip.position = "left") + 
  ggforce::facet_col(vars(cell), scales = 'free_y', strip.position = "left") + 
  scale_x_continuous(expand = c(0.001, 0.001)) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_blank(), 
        panel.spacing = unit(-.1, "line"), 
        strip.text = element_text(size = 12), 
        strip.background = element_rect(colour = NA, fill = NA), 
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0)) -> p10
p10
```

```{r fig.width=14, fig.height=9}
p10 <- p10 + scale_x_continuous(expand = c(0, 0))
p10 <- p10 + theme(plot.margin = unit(c(0, 0.5, 1, 1.9), "line"))
pe <- cowplot::plot_grid(p6, p7, p10, p8, ncol = 1, rel_heights = c(6, 1.75, 2.5, 2.5))
pe
```

```{r fig.width=14, fig.height=6}
p10_2 <- p10 + scale_y_sqrt()
cowplot::plot_grid(p10_2, p8, ncol = 1, rel_heights = c(3, 1.75))
```


```{r fig.width=14, fig.height=9.2}
pe_2 <- cowplot::plot_grid(p6, p7, p10_2, p8, ncol = 1, rel_heights = c(6, 1.75, 2.5, 2.5))
pe_2

ggsave(pe_2, filename = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/candi/TC/Capg.tc_tablus.pdf", width = 10, height = 8)

```

# d
```{r}
cell_de3 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/cell.expre_tss.filter.Rds")

cell_de3$end <- cell_de3$start
cell_de3[cell_de3$Group == "HSC",] %>% makeGRangesFromDataFrame() -> hsc_gr
lapply(cellorder, function(i){
  cell_de3[cell_de3$Group == i,] %>% makeGRangesFromDataFrame() -> gr2
  dist_result <- distanceToNearest(gr2, hsc_gr)
  valid <- dist_result[abs(mcols(dist_result)[,1])<=50,]
  data.frame(type=i,
             "Conserve"=round(length(unique(queryHits(valid))) / length(gr2)*100,0),
             "Loss"=round(length(setdiff(seq_along(hsc_gr),subjectHits(valid)))/ length(gr2)*100,0),
             "Gain"=round(length(setdiff(seq_along(gr2),queryHits(valid)))/length(gr2)*100,0)) 
}) %>% do.call(rbind,.) ->  plot


plot[plot$type=="HSC","Loss"] <- 0

cell_de3[cell_de3$Group == "HSC",]$Gene -> ctl
lapply(cellorder, function(i){
  cell_de3[cell_de3$Group == i,]$Gene -> ctl2
  data.frame(type=i,
             "Conserve"= length(intersect(ctl,ctl2)),
             "Loss"=length(setdiff(ctl,ctl2)),
             "Gain"=length(setdiff(ctl2,ctl)))
}) %>% do.call(rbind,.) -> plot
```
```{r fig.height=4, fig.width=8}
plot$Loss <- -plot$Loss
long_data <- plot %>% reshape2::melt()

level <- c("Loss","Conserve","Gain")
col <- c("#b7d0d7","#eaa08e","#df7359")
names(col) <- level
```


```{r fig.height=4, fig.width=8}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/0615_gain_loss_all_tss.pdf", height = 4, width = 5)

ggplot(long_data,
       aes(x=factor(type, levels=(cellorder)),
           y=value,
           fill=factor(variable, levels = rev(level))))+
geom_bar(stat = "identity",
         position ="stack",width = 0.9 )+
  theme_classic()+
  scale_fill_manual(values = col, name='')+
  ylab("Percentage of TSSs (%)")+ggtitle("Gain and Loss of TSS events from HSC")+
  scale_y_continuous()+
  xlab(NULL)+
 theme(legend.position = "top",
       legend.text = element_text(size = 12),
       axis.text.y = element_text(size = 12),
       axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
       axis.title  = element_text(size = 14)) +
  guides(fill = guide_legend( nrow = 1)) 


dev.off()
```

# e
```{r}
obj_df_new <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/tssmulti_one_proxiaml.Rds")
obj_df_new$type1 <- unlist(stringr::str_split_fixed(obj_df_new$type,"@",2))[,1]
obj_df_new$type2 <- unlist(stringr::str_split_fixed(obj_df_new$type,"@",2))[,2]
obj_df_new[obj_df_new$type1=="exist",] -> obj2


obj2$type2 <- ifelse(obj2$type2=="0","Distal",ifelse(obj2$type2=="1","Proximal",obj2$type2))
obj2 <- obj2[obj2$type2 !=0.5,]
```
```{r fig.height=3, fig.width=2}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/lineage_multi_one2.pdf", height = 3, width = 2)
table(obj2[obj2$start1!=obj2$start2,]$type2) %>% data.frame() ->res
res$Var1 <- factor(res$Var1, levels =  c("Proximal","Distal"))
ggplot(res,aes(x=Var1,y=Freq, fill = Var1))+
  geom_bar(stat = "identity",position =position_dodge(0.7),width = 0.7)+
  theme_classic()+
  labs(x="",y="Number of genes",title = "Lost ATSs \nfrom >1 to 1")+
  scale_fill_brewer(palette = "Pastel1")+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 14,hjust = .5))
dev.off()
```

# f
```{r}
library(ChIPseeker)
library(org.Mm.eg.db)

txdb <- makeTxDbFromGFF("/mnt/raid66/Personal_data/xuzijie/ref/MusMus/release101/Mus_musculus.GRCm38.101.chr.gtf", format = "gtf")
obj_gr_lst <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/phastcon_lost_gr.Rds")

lapply(obj_gr_lst, function(peaks){
  peak_anno <- annotatePeak(
    peaks,
    addFlankGeneInfo = T,
    genomicAnnotationPriority = c("Promoter", "5UTR","Exon", "3UTR", "Intron","Downstream", "Intergenic"),
    tssRegion = c(-50, 50),  # 定义启动子区域范围
    TxDb = txdb,  # 基因组注释数据库
    annoDb = "org.Mm.eg.db",      # 用于基因ID转换
    ignoreDownstream=F
    )
  peak_anno
}) -> peak_info

```

```{r}
col <- c("#a3c8dc","#92c875","#f19232","#a95728")
names(col) <- c("Promoter","5' UTR","Exon","Intron")
data.frame(class=c(rep("Proximal",4),rep("Distal",4)),
                            type=c(rep(names(col), 2)),
                            per=c(60.4545455,10.4545455, 29.09091,0,
                                  93.4782609,4.3478261,1.4492754,0.7246377), stringsAsFactors = F) -> plot_new

```

```{r fig.height=4, fig.width=4}
plot_new$type <- factor(plot_new$type, levels = names(col)) 
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/1202_lost_regions.pdf", height = 4, width = 4)

  ggplot(plot_new,
         aes(x=factor(class, levels=c("Proximal","Distal")),
             y=per,
             fill=reorder(type,-per),
             label=per))+
  geom_bar(stat = "identity",
           position ="fill",alpha=0.9,width = 0.8 )+
    theme_classic()+
    scale_fill_manual(values = c(col))+
    ylab("Percentage of TSSs (%)")+
    scale_y_continuous()+#百分比加这个
    xlab(NULL)+
   theme(legend.position = "right",
          legend.title = element_blank(),
         axis.title = element_text(size = 14),
         legend.text = element_text(size = 12),
         axis.text = element_text(size = 12))+
    guides(fill = guide_legend( nrow = 6))
dev.off()

```


# g
```{r}
obj_df_new <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/tssmulti_one_proxiaml.Rds")
obj_df_new$type1 <- unlist(stringr::str_split_fixed(obj_df_new$type,"@",2))[,1]
obj_df_new$type2 <- unlist(stringr::str_split_fixed(obj_df_new$type,"@",2))[,2]
obj_df_new[obj_df_new$type1=="exist",] -> obj2


obj2$type2 <- ifelse(obj2$type2=="0","Distal",ifelse(obj2$type2=="1","Proximal",obj2$type2))
obj2 <- obj2[obj2$type2 !=0.5,]
obj3 <- obj2[obj2$start1 != obj2$start2,]
obj3$chr <- mapply(function(y) y[1], strsplit(mapply(function(x) x[2], strsplit(obj3$TSS1, "@",2)), ":",2))
obj3$strand <- mapply(function(y) y[3], strsplit(mapply(function(x) x[2], strsplit(obj3$TSS1, "@",2)), ":",3))
obj3[,c("TSS1","start1","type2","chr","strand")] -> obj4
colnames(obj4)[2] <- "start"
obj4$end <- obj4$start
makeGRangesFromDataFrame(obj4, keep.extra.columns = T) -> obj_gr
split(obj_gr, obj_gr$type2) -> obj_gr_lst

lapply(c(1:2),function(i){
  print(i)
  phast <- get_phast(g_rds = obj_gr_lst[[i]],
                    bw_path = "/mnt/raid66/Personal_data/xuzijie/ref/MusMus/UCSC/mm10.60way.phastCons.bw",
                    seq = 250,
                    label = names(obj_gr_lst)[i])
  phast
})%>% do.call(rbind,.)->utr_lst


utr_lst[,  .( mean2 = mean(y)), by=ss]
utr_lst[,  .( mean2 = median(y)), by=ss]
```

```{r fig.width = 3.5,fig.height=3}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/phastcon_lost_line.pdf",height = 3,width = 3.5)

ggplot(utr_lst, aes(x = Pos, y = y,color = factor(ss,levels =  c("Proximal","Distal")))) + 
    geom_line(size = 0.3) +
  scale_color_manual(values = c("#f9a4a0","#8dbfe0"))+
  ylab("Average phastCons score")+xlab("Distance to TSS")+
    scale_x_continuous(breaks = seq(-250,250,125))+
    theme_classic()+
    theme(axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 8),
          axis.title =  element_text(size = 14),
          legend.text=element_text(size=14),
          plot.title = element_text (hjust = 0.5),
          legend.title = element_blank(),
          strip.text = element_text(size = 14),
          legend.position = c(0.8,0.7))

dev.off()
```

# h


# j

```{r}
library(clusterProfiler)
library(org.Mm.eg.db)
gain$type <- "Gained"
gain$ID <- mapIds(org.Mm.eg.db,
                  keys = gain$Group.2,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals="first")
ora_res <- enrichGO(gene = gain$ID,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENTREZID",#这里指定ID类型
                   ont = "BP", # "BP", "MF", "CC" 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   minGSSize = 5,# 最少的基因数量
                   maxGSSize = 300, # 最大的基因数量
                   readable = T # 把ENTREZID转换为SYMBOL
                   ) 
formu_res=clusterProfiler::simplify(ora_res,
                                         cutoff=0.5,
                                         by="p.adjust",
                                         select_fun=min)

saveRDS(formu_res, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/lineage.go.ats_gain.Rds")
```


```{r}
lost$type <- "Lost"
lost$ID <- mapIds(org.Mm.eg.db,
                  keys = lost$Group.2,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals="first")
ora_res <- enrichGO(gene = lost$ID,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENTREZID",#这里指定ID类型
                   ont = "BP", # "BP", "MF", "CC" 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   minGSSize = 5,# 最少的基因数量
                   maxGSSize = 300, # 最大的基因数量
                   readable = T # 把ENTREZID转换为SYMBOL
                   ) 
formu_res=clusterProfiler::simplify(ora_res,
                                         cutoff=0.5,
                                         by="p.adjust",
                                         select_fun=min)

saveRDS(formu_res, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/lineage.go.ats_lost.Rds")
```

```{r fig.height=4, fig.width=5}
ats_go_lost <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/lineage.go.ats_lost.Rds")
ats_go_gain <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/lineage.go.ats_gain.Rds")

rbind(ats_go_lost,ats_go_gain) -> res
res$labelx=rep(0,nrow(res))
res$labely=seq(nrow(res),1)
res$type <- factor(res$type, levels = (c("Lost","Gain")))
```


```{r fig.height=4, fig.width=5}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/0609_lineage.go.ats_lost.pdf", height = 4, width = 6)
ggplot(res, aes(x=factor(Description,levels =  (c(res$Description))), 
                               y=ifelse(type=="Gain",log10(pvalue),-log10(pvalue)), fill=type))+
  scale_fill_manual(values = (c("Lost"="#62c0dd", "Gain"="#d184aa")))+
  coord_flip()+
  geom_bar(stat="identity",
           alpha=0.5,
           width = 0.8) + 
  geom_text(data = subset(res,type=="Gain"),aes(x=labely,
                y=labelx,
                label = Description),
            size=3.5,hjust=1)+
  geom_text(data = subset(res,type=="Lost"),
            aes(x=labely,
                y=labelx,
                label = Description),
            size=3.5,hjust=0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none")+
  labs(x="", y="-Log10 (P value)", title = "GO (BP) of TSS host genes")

dev.off()

```

# k

```{r}
seu <- readRDS("/mnt/raid61/Personal_data/xuzijie/task/07ATS/02result/thesis/mmBM/V4/01quant/mmBM.DF.ATS.Rds")
seu@assays$RNA@counts -> df
df[rowSums(df)>0,] -> df2
write.csv(t(as.matrix(df2)),file = "/mnt/raid61/Personal_data/xuzijie/task/07ATS/02result/thesis/mmBM/V4/07SCENIC/sce_exp.csv")
```


```{bash}
cd /mnt/raid61/Personal_data/xuzijie/task/07ATS/02result/thesis/mmBM/V4/07SCENIC
python trans.py

cd /mnt/raid61/Personal_data/xuzijie/task/07ATS/02result/thesis/mmBM/V4/07SCENIC
pyscenic grn --num_workers 15 \
--sparse \
--method grnboost2 \
--output sce.adj.csv \
sce.loom \
/mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/SCENIC/mm_mgi_tfs.txt # tfs names, 

pyscenic ctx --num_workers 15 \
--output sce.regulons.csv \
--expression_mtx_fname sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/SCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
sce.adj.csv \
/mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/SCENICmm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

pyscenic aucell --num_workers 40 \
--output sce_SCENIC_xzj.loom \
sce.loom \
sce.regulons.csv

```

```{r fig.height=2.5, fig.width=3}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/06TF/gene_exp_dotplot.pdf", height = 2.5, width = 3)

DotPlot(seu, features = rev(c("Rest","Yy1","Spi1","Zmiz1","Stat3")), group.by = "Lineages")+coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x='',y='')
dev.off()
```

# l
```{r}
library(SCENIC)
meta <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V3/00quant/DF_UTR.colData.Rds") %>% data.frame(.)

cellanno_lineage <- meta %>%dplyr::select(group)

sce_SCENIC <- open_loom("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/covid19/sce_SCENIC2.loom")

regulons.500bp.tss <- get_regulons(sce_SCENIC, column.attr.name="Regulons") # 行为regulon，列为基因的ont-hot矩阵，行的加和即为regulonsToGeneLists的结果

regulonsGene.500bp.tss <- regulonsToGeneLists(regulons.500bp.tss)#每个regulon靶向的基因
regulonAUC.500bp.tss <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC') # 每个cellbc的regulon活性

regulonAUC.500bp.tss2 <- regulonAUC.500bp.tss[,rownames(cellanno_lineage)]

RSS.500bp.tss.lineage <- SCENIC::calcRSS(AUC = AUCell::getAUC(regulonAUC.500bp.tss2),cellAnnotation = cellanno_lineage[colnames(AUCell::getAUC(regulonAUC.500bp.tss2)),"group"])

```

```{r}

rssPlot <- SCENIC::plotRSS(
  RSS.500bp.tss.lineage,
  zThreshold = 0,
  cluster_columns = FALSE,
  order_rows = TRUE,
  varName = "group"
)
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/06TF/RSS_dotplot.pdf", height = 2.5, width = 2.5)

rssPlot$plot +
  theme(axis.text.x = element_text(angle = 45))
dev.off()

```


# m
```{r}
library(dplyr)
library(tidyr)

proj <- readRDS("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/mouse/prog_3TF_TSS.Rds")
proj$id <- ifelse(proj$id=="Progenitor","HSPCs","Mature")
proj$id <- factor(proj$id, levels = lineage)
mature <- readRDS("/mnt/raid66/Personal_data/wanghuiying/04.ATS/tss_motif_analyse/scenic/mouse/mature_2TF_TSS.Rds")
mature$id <- ifelse(mature$id=="Progenitor","HSPCs","Mature")
mature$id <- factor(mature$id, levels = lineage)

```

```{r fig.height=3, fig.width=4.5}
library(ggtext)

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/06TF/mature_tf_dotplot.pdf", height = 3, width = 4.5)
ggplot(mature,
       aes(x = factor(id, levels = rev(lineage)),
           y = features.plot)) +
  geom_point(
    aes(size = pct.exp, color = avg.exp),
    alpha = 0.8
  ) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue",
    name = "Average \nExpression"
  ) +
  scale_size(range = c(1, 6), name = "Percentage \nExpression") +
  new_scale_color() +
  geom_point(
    data = label_df,
    aes(
      x = 1,
      y = features.plot,
      color = TF
    ),
    alpha = 0,
    size = 3
  ) +
  scale_color_manual(
    values = tf_colors,
    name = "Mature-dominant TF"
  ) +
  scale_y_discrete(labels = label_p) +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1,size = 10),
    axis.text.y = element_text(size=10),
    legend.position = "right"
  ) +
  labs(x = NULL, y = "")+
  guides(
  color = guide_legend(
    override.aes = list(alpha = .6, size = 4)
  )
)+coord_flip()
dev.off()

```

```{r fig.height=8, fig.width=4}
proj$TF <- factor(
  proj$TF,
  levels = c("Rest", "Spi1", "Yy1")
)

tf_colors <- c(
  "Rest" = "#D55E00",
  "Spi1"  = "#7570B3",
  "Yy1"   = "#1b9e77"
)

label_df <- proj %>%
  distinct(features.plot, TF) %>%
  mutate(
    label_colored = paste0(
      "<span style='color:",
      tf_colors[TF],
      "'>",
      features.plot,
      "</span>"
    )
  )

label_p <- label_df$label_colored
names(label_p) <- label_df$features.plot
scales::show_col(brewer.pal(4,"Dark2"))
```


```{r fig.height=3, fig.width=4.5}
library(ggtext)
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/06TF/HSPC_tf_dotplot.pdf", height = 3, width = 4.5)

ggplot(proj,
       aes(x = factor(id, levels = rev(lineage)),
           y = features.plot)) +
  geom_point(
    aes(size = pct.exp, color = avg.exp),
    alpha = 0.95
  ) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue",
    name = "Average \nExpression"
  ) +
  scale_size(range = c(1, 6), name = "Percentage \nExpression") +
  new_scale_color() +
  geom_point(
    data = label_df,
    aes(
      x = 1,
      y = features.plot,
      color = TF
    ),
    alpha = 0,
    size = 3
  ) +
  scale_color_manual(
    values = tf_colors,
    name = "HSPCs-dominant TF"
  ) +
  scale_y_discrete(labels = label_p) +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1,size = 10),
    axis.text.y = element_text(size=10),
    # legend.text = element_markdown(size = 10),
    legend.position = "right"
  ) +
  labs(x = NULL, y = "")+
  guides(
  color = guide_legend(
    override.aes = list(alpha = .6, size = 4)
  )
)+coord_flip()
dev.off()



```

