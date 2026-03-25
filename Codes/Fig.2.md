---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
 
# load
```{r warning=F}

suppressMessages({
  source("/mnt/raid66/Personal_data/xuzijie/task/07ATS/01script/thesis/00functions.R")
  library(scATS)
  library(ggplot2)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(magrittr)
  library(RColorBrewer)
  library(parallel)
})

c("#FC8D62","#66C2A5","#8DA0CB") -> col
names(col) <- c("scATS","SCAFE","CamoTSS")

c("#FC8D62","#66C2A5","#8DA0CB","#E78AC3" ,"#a9b2b6") -> col
names(col) <- c("scATS","SCAFE","CamoTSS","CAGE","Random")

c("#FC8D62","#66C2A5","#8DA0CB", "#ecd276","#da954f") -> col
names(col) <- c("scATS","SCAFE","CamoTSS","TSSr","scTSS")

c("#FC8D62","#66C2A5","#8DA0CB","#E78AC3" ,"#ED646B","#a9b2b6" ,"#ecd276","#da954f") -> col
names(col) <- c("scATS","SCAFE","CamoTSS","GENCODE","CAGE","Random","TSSr","scTSS")



```

# a-b
```{r fig.height=3, fig.width=4}
tss_df2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATS_obj_5_0105.Rds")

len <- lapply(tss_df2, length) %>% data.frame(.)  %>% t() %>% data.frame(num=.) %>% tibble::rownames_to_column(.,var="type")
len[1,]$num <- 10195
ggplot(len, aes(x = factor(type, levels = rev(c("scATS","CamoTSS","scTSS","SCAFE","TSSr"))), y = num,fill=type)) +
  geom_bar(stat = "identity",position = "dodge", color="black", alpha=.5) +
  theme_classic() +
  labs(
    x = "",
    y = "Number of TSSs",
    fill = "",
    title = ""
  ) +
  geom_text(aes(label=num),
            size=4.5,
            position = position_dodge())+
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 12),
    legend.position = "top"
  ) +coord_flip()+
  scale_fill_manual(values = col)->p

```

```{r}
tss_df2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATS_obj_5_0105.Rds")

get_venn_id <- function(df){
  gr_extended <- resize(df, width = width(df)+50, fix = "center")
  gr_extended
}

scats <- get_venn_id(tss_df2$scATS)#B
scafe <- get_venn_id(tss_df2$SCAFE)#B
camotss <- get_venn_id(tss_df2$CamoTSS)#C
tssr <- get_venn_id(tss_df2$TSSr)#B
sctss <- get_venn_id(tss_df2$scTSS)#C


other <- list(scats, scafe, camotss,tssr,sctss)
names(other) <- c("scATS", "SCAFE", "CamoTSS","TSSr","scTSS")
len <- lapply(other, length)

lapply(1:5, function(i){
    
  anchor <- other[[i]]
  set <- other
  set[[i]] <- NULL
  detect_mat <- sapply(
    set,
    function(s) {
      ov <- findOverlaps(anchor, s)
      detected <- rep(0, length(anchor))
      detected[queryHits(ov)] <- 1
      detected
    }
  )

  overlap_count <- rowSums(detect_mat)

  df_new <- data.frame(
  group = c("0", "1", "2", "3","4"),
    Freq = c(
      sum(overlap_count == 0),
      sum(overlap_count == 1),
      sum(overlap_count == 2),
      sum(overlap_count == 3),
      sum(overlap_count == 4)
    )
  )


  df_new$type <- names(other)[i]
  df_new$len <- len[[names(other)[i]]]
  df_new$prop <- df_new$Freq / df_new$len

  df_new
}) %>% do.call(rbind, .) -> df_lst
df_lst
```

```{r fig.height=3, fig.width=4}
color <- c("grey","#faddbc","#da944a","#733915","#340E05")

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATS_obj_5_0105_validate.pdf", height = 3, width = 4)
ggplot(df_lst, aes(x = factor(type, levels = rev(c("scATS","CamoTSS","scTSS","SCAFE","TSSr"))), y = prop, fill = factor(group, levels = (c(0:4))))) +
  geom_bar(stat = "identity",position = "fill", color="black", alpha=.4) +
  theme_classic() +
  labs(
    x = "",
    y = "TSSs validated by tools (%)",
    fill = "Identified by",
    title = ""
  ) +
  geom_text(aes(label=round(prop*100,0)),
            size=4.5,
            position = position_fill(vjust=0.5))+
  theme(
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 12),
    legend.position = "top"
  ) +coord_flip()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,"Blues"), name="")->p4

```


```{r}

tss_df2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATS_obj_5_0105.Rds")

gtf_dis_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/distance/all5_gtf.Rds")
reftss_dis_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/distance/all5_reftss.Rds")


gtf_dis_lst[[1]]$label <- paste0(gtf_dis_lst[[1]]$gene_name,"@",gtf_dis_lst[[1]]$TSS)
gtf_dis_lst[[3]]$label <- paste0("CamoTSS@",gtf_dis_lst[[3]]$seqnames,"@",gtf_dis_lst[[3]]$start,"@",gtf_dis_lst[[3]]$strand)
gtf_dis_lst[[2]]$label <- paste0("SCAFE@",gtf_dis_lst[[2]]$seqnames,"@",gtf_dis_lst[[2]]$start,"@",gtf_dis_lst[[2]]$width,"@",gtf_dis_lst[[2]]$strand)
gtf_dis_lst[[4]]$label <- paste0("TSSr@",gtf_dis_lst[[4]]$seqnames,"@",gtf_dis_lst[[4]]$start,"@",gtf_dis_lst[[4]]$width,"@",gtf_dis_lst[[4]]$strand)
gtf_dis_lst[[5]]$label <- paste0("scTSS@",gtf_dis_lst[[5]]$seqnames,"@",gtf_dis_lst[[5]]$start,"@",gtf_dis_lst[[5]]$width,"@",gtf_dis_lst[[5]]$strand)



reftss_dis_lst[[1]]$label <- paste0(reftss_dis_lst[[1]]$gene_name,"@",reftss_dis_lst[[1]]$TSS)
reftss_dis_lst[[3]]$label <- paste0("CamoTSS@",reftss_dis_lst[[3]]$seqnames,"@",reftss_dis_lst[[3]]$start,"@",reftss_dis_lst[[3]]$strand)
reftss_dis_lst[[2]]$label <- paste0("SCAFE@",reftss_dis_lst[[2]]$seqnames,"@",reftss_dis_lst[[2]]$start,"@",reftss_dis_lst[[2]]$width,"@",reftss_dis_lst[[2]]$strand)
reftss_dis_lst[[4]]$label <- paste0("TSSr@",reftss_dis_lst[[4]]$seqnames,"@",reftss_dis_lst[[4]]$start,"@",reftss_dis_lst[[4]]$width,"@",reftss_dis_lst[[4]]$strand)
reftss_dis_lst[[5]]$label <- paste0("scTSS@",reftss_dis_lst[[5]]$seqnames,"@",reftss_dis_lst[[5]]$start,"@",reftss_dis_lst[[5]]$width,"@",reftss_dis_lst[[5]]$strand)


lapply(1:5, function(i){
  anno_gtf <- gtf_dis_lst[[i]][abs(gtf_dis_lst[[i]]$distance ) <= 50,]$label
  anno_ref <- reftss_dis_lst[[i]][abs(reftss_dis_lst[[i]]$distance ) <= 50,]$label
  both <- intersect(anno_gtf, anno_ref)
  only_gtf <- setdiff(anno_gtf, both)
  only_ref <- setdiff(anno_ref, both)
  novel <- setdiff(gtf_dis_lst[[i]]$label,c(anno_gtf, anno_ref))
  data.frame(type = names(tss_df2)[i],name=c("Novel","GENCODE and refTSS","GENCODE or refTSS"),
             num = c(length(novel), length(both), c(length(only_ref)+ length(only_gtf)))) -> plot
  plot$perc <- round(plot$num *100 /sum(plot$num),0)
  plot
}) %>% do.call(rbind,.) -> plot
```


```{r fig.height=3, fig.width=4}
color <- c("grey","#faddbc","#da944a","#733915","#340E05")

ggplot(plot,
       aes(x=factor(type, levels = rev(c("scATS","CamoTSS","scTSS","SCAFE","TSSr"))),
           y=num,
           fill=factor(name, levels = rev(c("GENCODE and refTSS","GENCODE or refTSS","Novel"))),
           label=perc))+
geom_bar(stat = "identity",
         position ="fill",color="black",alpha=0.5 )+
  theme_classic()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,"Blues")[1:3], name="")+
  ylab("TSSs validated by database (%)")+
  scale_y_continuous()+#百分比加这个
  xlab(NULL)+
  geom_text(aes(label=perc),
            size=4.5,
            position = position_fill(vjust=0.5))+ #调整标签位置
 theme(legend.position = "top",
       legend.text = element_text(size = 12),
       axis.text = element_text(size = 12),
       axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
       axis.title.x  = element_text(size = 12)) +
  coord_flip()->p2 

```


```{r}
TNS_TSS <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/TNG/all_TSS_v2.Rds")
lapply(TNS_TSS, function(df){
  data.frame(df, stringsAsFactors = F) -> df
  df$seqnames <- as.character(df$seqnames)
  df$seqnames <- ifelse(df$seqnames %in% c(as.character(c(1:21)),"MT","M","X","Y"), paste0("chr", df$seqnames), df$seqnames)
  df %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
}) -> TNS_TSS
names(TNS_TSS) <- c("ONT", "PacBio", "scONT","Smartseq")


tss_df2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATS_obj_5_0105.Rds")

get_venn_id <- function(df){
  gr_extended <- resize(df, width = width(df)+50, fix = "center")
  gr_extended
}


ont <- get_venn_id(TNS_TSS$ONT) # A
pacbio <- get_venn_id(TNS_TSS$PacBio) # A
scONT <- get_venn_id(TNS_TSS$scONT) # A
Smartseq <- get_venn_id(TNS_TSS$Smartseq) # A


scats <- get_venn_id(tss_df2$scATS)#B
scafe <- get_venn_id(tss_df2$SCAFE)#B
camotss <- get_venn_id(tss_df2$CamoTSS)#C
tssr <- get_venn_id(tss_df2$TSSr)#B
sctss <- get_venn_id(tss_df2$scTSS)#C

```


```{r}
sets <- list(
  CamoTSS = camotss,
  ONT = ont,
  PacBio = pacbio
)

other <- list(scats, scafe, camotss,tssr,sctss)
names(other) <- c("scATS", "SCAFE", "CamoTSS","TSSr","scTSS")
len <- lapply(other, length)

lapply(1:5, function(i){
    
  anchor <- other[[i]]

  detect_mat <- sapply(
    list(ONT=ont, PacBio=pacbio, scONT=scONT),
    function(s) {
      ov <- findOverlaps(anchor, s)
      detected <- rep(0, length(anchor))
      detected[queryHits(ov)] <- 1
      detected
    }
  )

  overlap_count <- rowSums(detect_mat)

  df_new <- data.frame(
  group = c("0", "1", "2", "3"),
    Freq = c(
      sum(overlap_count == 0),
      sum(overlap_count == 1),
      sum(overlap_count == 2),
      sum(overlap_count == 3)
    )
  )


  df_new$type <- names(other)[i]
  df_new$len <- len[[names(other)[i]]]
  df_new$prop <- df_new$Freq / df_new$len

  df_new
}) %>% do.call(rbind, .) -> df_lst

```

```{r fig.height=3, fig.width=4}
color <- c("grey","#faddbc","#da944a","#733915")
levels <- c("scATS","SCAFE","CamoTSS","TSSr","scTSS")
# color <- pal_npg("nrc", alpha = 0.7)(9)
names(color) <- c(0:3)

ggplot(df_lst, aes(x = factor(type, levels = rev(c("scATS","CamoTSS","scTSS","SCAFE","TSSr"))), y = prop, fill = factor(group, levels = (c(0:3))))) +
  geom_bar(stat = "identity",position = "fill", color="black", alpha=.5) +
  theme_classic() +
  labs(
    x = "",
    y = "TSSs validated by long-read datasets (%)",
    fill = "Identified by",
    title = ""
  ) +
  theme(
    axis.text = element_text(size = 12),
    text = element_text(size = 12),
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) +coord_flip()+
   scale_fill_manual(values = RColorBrewer::brewer.pal(5,"Blues")[1:4], name="")+
  geom_text(aes(label=round(prop*100,0)),
            size=4.5,
            position = position_fill(vjust=0.5)) -> p3

```

```{r fig.height=3, fig.width=10}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/barplot_lst_2.pdf", height = 3, width = 10)

cowplot::plot_grid(p,p4,p2,p3, nrow = 1, rel_widths = c(1.1,0.9,0.65,0.8), rel_heights = c(1,0.2,1,1))
dev.off()
```

# c
```{r fig.width=9,fig.height=3}
library(dplyr)
library(openxlsx)
tss_df2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/distance/all5_gtf.Rds")
names(tss_df2) <- names(col)
ref <- readRDS("/mnt/raid61/Personal_data/xuzijie/task/07ATS/02result/tangchao/totalseq/02scATS_quant/ref.Rds")
```


```{r}
breaks=10;max=100;min=0
df_tc <- pre_fun(tss_df2$scATS,ref$gtfTSS[ref$gtfTSS$site %in% tss_df2$scATS$site,], exten = seq(min,max,breaks)) %>% mutate(type="scATS")
tss_df2$SCAFE$TSS <-paste0("SCAFE@",tss_df2$SCAFE$seqnames,"@",tss_df2$SCAFE$start,"@",tss_df2$SCAFE$width,"@",tss_df2$SCAFE$strand)
df_scafe <- pre_fun(tss_df2$SCAFE,ref$gtfTSS[ref$gtfTSS$site %in% tss_df2$SCAFE$site,], exten = seq(min,max,breaks)) %>% mutate(type="SCAFE")

tss_df2$CamoTSS$TSS <-paste0("CamoTSS@",tss_df2$CamoTSS$seqnames,"@",tss_df2$CamoTSS$start,"@",tss_df2$CamoTSS$width,"@",tss_df2$CamoTSS$strand)
df_camo <- pre_fun(tss_df2$CamoTSS,ref$gtfTSS[ref$gtfTSS$site %in% tss_df2$CamoTSS$site,], exten = seq(min,max,breaks)) %>% mutate(type="CamoTSS")

tss_df2$TSSr$TSS <-paste0("TSSr@",tss_df2$TSSr$seqnames,"@",tss_df2$TSSr$start,"@",tss_df2$TSSr$width,"@",tss_df2$TSSr$strand)
df_tssr <- pre_fun(tss_df2$TSSr,ref$gtfTSS[ref$gtfTSS$site %in% tss_df2$TSSr$site,], exten = seq(min,max,breaks)) %>% mutate(type="TSSr")

tss_df2$scTSS$TSS <-paste0("scTSS@",tss_df2$scTSS$seqnames,"@",tss_df2$scTSS$start,"@",tss_df2$scTSS$width,"@",tss_df2$scTSS$strand)
df_sctss <- pre_fun(tss_df2$scTSS,ref$gtfTSS[ref$gtfTSS$site %in% tss_df2$scTSS$site,], exten = seq(min,max,breaks)) %>% mutate(type="scTSS")

levels <- rev(names(col))

df <- rbind(
  df_tc,df_scafe,df_camo, df_tssr, df_sctss)
df$type <- factor(df$type,levels = rev(levels))
df$class <- factor(df$class,levels = c("Precision","Recall","Fscore"))


```


```{r fig.width=7,fig.height=2.8}
c("#FC8D62","#66C2A5","#8DA0CB", "#ecd276","#da954f") -> col
names(col) <- c("scATS","SCAFE","CamoTSS","TSSr","scTSS")

pdf(paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/","v8_gtf_pre_filter_1206.pdf"),width = 7,height = 2.8)
my_plot(df, title = NULL) +facet_wrap(~ class)+
  scale_color_manual(values = (col))+
  theme(
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        legend.position = c("top"),
        strip.background =  element_rect(colour="black",
                                        fill="white"),
          strip.text = element_text(size = 12, color = "black")
        )+
  scale_x_continuous(breaks = seq(0,100,20))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  labs(y="The ratio of matched TSSs", x = "Cutoff (bp)")

dev.off()

```

# e
```{r echo=FALSE, message=FALSE}
suppressMessages({
  library(magrittr)
  library(data.table)
  library(ggplot2)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  # for GSEA analysis
  library(msigdbr)
  
  # For motif enrichment analysis
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  # MOFA
  library(MOFA2)
  library(MOFAdata)
  library(RColorBrewer)
})

cellorder <- c("HSC" ,"CMP", "GMP", "EP", "MKP","BasoP", "Baso", "Neu", "T cell", "B cell", "NK")
c(brewer.pal(11, "Paired")) -> color
color[11] <- "#b3b36b"
names(color) <- cellorder

col <- c("Pre-mature", "Mature")
color2 <- c('#F4BC58','#8DA0CB')
names(col) <- color2
data.frame(col)
```


```{r}
seu <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/bm.DF.ATS.Rds")

seu@meta.data %>% tibble::rownames_to_column(., var = "sample") -> bm_meta

MOFAobject0 <- create_mofa(seu, assays = c("RNA", "Protein", "ats_Theta"))

data_opts <- get_default_data_options(MOFAobject0)

model_opts <- get_default_model_options(MOFAobject0)
model_opts$num_factors <- 15

train_opts <- get_default_training_options(MOFAobject0)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42


MOFAobject0 <- prepare_mofa(MOFAobject0,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
plot_data_overview(MOFAobject0) 
MOFAobject <- run_mofa(MOFAobject0, outfile="/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/MOFA2_bm_theta_0424.hdf5")
saveRDS(MOFAobject,"/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/MOFA2_bm_theta_0424.Rds")

```

```{r fig.height=3, fig.width=3}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/bm.factors.decomposition2.1.pdf",height = 3, width = 2)
plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]]+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))
dev.off()
```

# f

```{r fig.height=3, fig.width=4}

factors2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/MOFA2_bm_theta_0424_factor.Rds")
res3 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/multi.ats.num.Rds")

boxplot_data <- factors2[factors2$factor=="Factor2",]
line_data <- res3
colnames(line_data)[1] <- "wsnn_cell_type"

boxplot_data$wsnn_cell_type <- factor(boxplot_data$wsnn_cell_type, levels = cellorder)
line_data$wsnn_cell_type <- factor(line_data$wsnn_cell_type, levels = cellorder)
boxplot_data$Lineages <- ifelse(boxplot_data$Lineages=="Mature","Mature","HSPCs")

scaling_factor <- max(boxplot_data$value) / max(line_data$x)
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/factor2_3.pdf", height = 3, width = 4)

ggplot(boxplot_data, aes(x = wsnn_cell_type)) + 
  geom_boxplot(aes(y = value,color =Lineages), alpha = 0.7) +
  geom_line(
    data = line_data, 
    aes(y = x * scaling_factor, group = 1), 
    color = "#e6855f", 
    linewidth = .5
  ) +
  geom_point(
    data = line_data,
    aes(y = x * scaling_factor),
    color = "#e6855f",
    size = 2
  ) +
  scale_y_continuous(
    name = "MOFA+ F2",
    sec.axis = sec_axis(~ . / scaling_factor, name = "Percentage of ATS (%)")
  ) +
  labs(
    title = "",
    x = ""
  ) +
  theme_test(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = .9, vjust = .9, size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        # axis.ticks.y.right = element_rect(color = "#e6855f"),
        legend.position = c(0.68,0.98))+
  geom_hline(yintercept=c(0), linetype="dotted", color = "#e6855f")+
  scale_color_manual(values = c("#da89b1","#bfb7da"), name = "")+
  guides(color = guide_legend( nrow = 1))
dev.off()
```

# g
```{r fig.height=3, fig.width=6}
library(Rtsne)

factors2 <- merge(MOFAobject@samples_metadata[,c("sample","wsnn_cell_type","Lineages")], factors, by="sample")
factors2$factor <- as.character(factors2$factor)
factors2[factors2$factor=="Factor2",] -> fact2
fact2$wsnn_cell_type <- as.character(fact2$wsnn_cell_type)

pseudo_bulk_df <- aggregate(fact2$value, by = list(cell_type=fact2$wsnn_cell_type), FUN = mean)
pseudo_bulk_df <- pseudo_bulk_df %>% tibble::column_to_rownames(.,var = "cell_type")

pca_result <- prcomp(pseudo_bulk_df, center = TRUE, scale. = F)

pca_df <- as.data.frame(pca_result$x)
pca_df$cell_type <- rownames(pca_df) # 将行名（细胞类型）添加为一列
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

```


```{r fig.height=4, fig.width=4}
num_samples <- nrow(pseudo_bulk_df)
# 计算最大可能的perplexity
max_perplexity <- floor((num_samples - 1) / 3)
# 设置一个安全值，但不超过一个合理的上限（例如10），并确保它至少为1
safe_perplexity <- max(1, min(10, max_perplexity - 1)) 


set.seed(123) # 设置随机种子以保证结果可重复
tsne_result <- Rtsne(
  as.matrix(pseudo_bulk_df),
  dims = 2,
  perplexity = safe_perplexity,
  verbose = FALSE,
  check_duplicates = FALSE # 对于低维度数据，可能存在重复行
)

# 将t-SNE结果转换为数据框
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c("tSNE_1", "tSNE_2")
tsne_df$cell_type <- rownames(pseudo_bulk_df) # 添加细胞类型信息
tsne_df$cell_type2 <- ifelse(tsne_df$cell_type %in% c("CMP", "BasophilP", "HSC", "GMP", "MKP", "EBP"), "HSPCs", "Mature")
```


```{r fig.height=3, fig.width=3}
col <- c("HSPCs", "Mature")
color2 <- c('#F4BC58','#8DA0CB')
names(col) <- color2

c(brewer.pal(11, "Paired")) -> color
color[11] <- "#b3b36b"
names(color) <- cellorder


pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/03mofa/tsne_factor2.pdf", height = 3,width = 3)
ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, color = wsnn_cell_type, label = wsnn_cell_type)) +
  geom_point(size = 5, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 4, max.overlaps = Inf) +
  labs(
    title =  "MOFA+ F2",
    x = "tSNE 1",
    y = "tSNE 2",
    color = "Cell Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text = element_text(size=14),
    legend.position = "none"
  )+
  scale_color_manual(values = color)+
  xlim(-175,185)
dev.off()
```

# h

```{r}

library(data.table)
library(ggplot2)
library(ggpubr)
library(ICAS)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(GenomicAlignments)
library(ggbio)

gtfFile <- file.path("/mnt/raid66/Personal_data/xuzijie/ref/MusMus/release101/Mus_musculus.GRCm38.101.chr.gtf")
gtf <- as.data.table(rtracklayer::readGFF(gtfFile))
txdb <- makeTxDbFromGFF(gtfFile)

unique(exons(txdb)[queryHits(findOverlaps(exons(txdb), as("chr6:72544391-72562983:+", "GRanges"), maxgap = 1))])
gr <- range(unique(exons(txdb)[queryHits(findOverlaps(exons(txdb), as("chr6:72544391-72555272:+", "GRanges"), maxgap = 1))]))
gr <- keepSeqlevels(gr, as.character(runValue(seqnames(gr))))
```

```{r fig.width=14, fig.height=5}
gr.txdb <- biovizBase::crunch(txdb, which = gr)
gr.txdb <- split(gr.txdb, gr.txdb$tx_name)
names(gr.txdb) <- plyr::mapvalues(names(gr.txdb), gtf[type == "transcript", transcript_id], gtf[type == "transcript", transcript_name], warn_missing = F)
sk <- mapply(function(x) {
  y <- subset(x, type == "exon")
  paste(as.character(sort(ranges(y))), collapse = ",")
}, gr.txdb)
sk <- split(names(sk), sk)
names(sk) <- mapply(function(x) paste(sort(x), collapse = ", "), sk)

gr.txdb <- lapply(sk, function(x) gr.txdb[[which(names(gr.txdb) == x[1])]])
gr.txdb <- GRangesList(gr.txdb)

p8 <- ggbio::autoplot(gr.txdb)@ggplot + scale_x_continuous(expand = c(0, 0), limits = c(start(gr), end(gr))) + ggbio::theme_clear()

p8 <- p8 + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  geom_vline(xintercept = c(72544428,72549434,72549458), color="red", linetype="dashed")
p8
```

```{r}
CellMeta <- fread("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/sashimi/cellbc_allcell.txt", header = F)
colnames(CellMeta) <- c("CB", "cell")
CellMeta <- CellMeta[CB %in% CellMeta[, .N, CB][N == 1, CB]]

gr <- GRanges(seqnames(gr), IRanges(min(layer_scales(p8)$x$range$range), max(layer_scales(p8)$x$range$range)), strand = strand(gr))
library(parallel)
minMapQuality <- 0
param0 = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isNotPassingQualityControls = FALSE), 
                                 tag = "CB",
                                 which = gr, 
                                 mapqFilter = minMapQuality)
bamfile <- file.path("/mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/bam/01bam/merge.bam")
bam <- GenomicAlignments::readGAlignments(file = bamfile, param = param0)
bam <- subset(bam, !is.na(CB))
```

```{r}
bam <- subset(bam, CB %in% gsub("^merge_", "", CellMeta[[1]]))
mcols(bam)$cell <- plyr::mapvalues(mcols(bam)$CB, gsub("^merge_", "", CellMeta[[1]]), CellMeta[[2]], warn_missing = FALSE)


CovList <- mclapply(split(bam, mcols(bam)$CB), function(i) {
  as.integer(coverage(i)[[seqlevels(gr)]][start(gr):end(gr)])
}, mc.cores = 40)

Cov <- as.data.frame(do.call(rbind, CovList))
row.names(Cov) <- names(CovList)
colnames(Cov) <- start(gr):end(gr)
Cov <- Cov[apply(Cov, 1, max) > 0, ]
Cov <- Cov[apply(Cov, 1, function(x) quantile(x, 0.99)) >= 1, ]
Cov <- Cov / apply(Cov, 1, function(x) quantile(x, 0.99))

melt_Cov <- data.table::melt.data.table(data.table::as.data.table(Cov, keep.rownames = "Cell"), variable.name = "Pos", value.name = "Coverage")
melt_Cov[, Pos := as.numeric(as.character(Pos))]
melt_Cov <- merge(melt_Cov, unique(CellMeta[, .(Cell = gsub("^merge_", "", CB), cell)]), by = "Cell")

melt_Cov[, cell := factor(cell, levels = cellorder)]
data.table::setkey(melt_Cov, cell)
Cell_od <- unique(melt_Cov[, .(Cell, cell)])[, Cell]
melt_Cov[, Cell := factor(Cell, levels = Cell_od)]
melt_Cov[, cell := factor(cell, levels = cellorder)]
C2U <- melt_Cov[, max(Coverage), Cell][V1 > 0, as.character(Cell)]
```


```{r}
anncolor <- list(Cell = colorRampPalette(rev(brewer.pal(n = 7, name = "Set1")))(11))
names(anncolor$Cell) <- cellorder

ggplot(data = melt_Cov, 
       mapping = aes(x = Pos, y = Cell, fill = Coverage)) + 
  geom_tile(show.legend = FALSE) + 
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

```{r fig.width=10, fig.height=8}
yls <- layer_scales(p8)$y$range$range
p8 <- p8 + scale_y_continuous(limits = c(min(yls), max(yls) + 0.5))
p8 <- p8 + theme(plot.margin = unit(c(0, 0.5, 1, 6.1), "line"))
pe <- cowplot::plot_grid(p6, p8, ncol = 1, rel_heights = c(6, 3.5))

pdf( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/candi/TC/Capg.tc_3.pdf", width = 10, height = 8)
pe
dev.off()
```


# j
```{r}
library(parallel)
library(stringr)
linea_psi_fil <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/cell.expre_tss.filter_lineage.Rds")
order <- c("Progenitor", "Mature")
color <- c('#F4BC58','#8DA0CB')
names(order) <- color


lapply(split(linea_psi_fil, linea_psi_fil$Group),function(df){
  mclapply(split(df, df$Gene), function(ge){
    data.frame(Var1 = unique(ge$Gene),
               Freq = length(unique(ge$TSS)),
               Group = unique(ge$Group))
  }, mc.cores = 20) %>% do.call(rbind,.) -> res
  res
})%>% do.call(rbind,.)->cell_ats

aggregate(cell_ats$Freq, list(cell_ats$Group,cell_ats$Var1), sum) -> cell_df

names(which(table(cell_df$Group.2) == 2)) -> loci

cell_df[cell_df$Group.2 %in% loci,] -> cell_df2 #河流图对象

my_merge <- function(df1, df2){
  merge(df1, df2, by =c("Group.2"))
}
cell_df2$x2 <- ifelse(cell_df2$x >= 2,">1",cell_df2$x)

lapply(split(cell_df2, cell_df2$Group.1), function(df){
  colnames(df)[4] <- unique(as.character(df$Group.1))
  df[,1] <- NULL
  df
}) %>% Reduce(my_merge,.) -> cell_df3
colnames(cell_df3)[5] <- "Progenitor"
cell_df3$pattern <- paste(cell_df3[,"Progenitor"], cell_df3[,"Mature"],sep = "_")

aggregate(cell_df3$pattern, list(cell_df3$Progenitor, cell_df3$Mature, cell_df3$pattern), length) ->cell_df4
colnames(cell_df4)[1:2] <- c("Progenitor", "Mature")
cell_df4$Progenitor <- factor(cell_df4$Progenitor,levels = c("1",">1"))
cell_df4$Mature <- factor(cell_df4$Mature,levels = c("1",">1"))

```
```{r fig.height=4,fig.width=4}
library(ggalluvial)
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/lineage.ggalluvial.alltss.pdf", height = 4, width = 4)

ggplot(cell_df4,
       aes(y = x,
           axis1 = Progenitor, axis2 = Mature)) +
  geom_alluvium(aes(fill = Progenitor),
                width = 0, reverse = FALSE) +
  guides(fill = FALSE) +theme_classic()+
  geom_stratum(width = 1/4, reverse = FALSE, color = "black")+
  geom_text(stat = "stratum", infer.label = TRUE, reverse = FALSE,size=5) +
  scale_x_discrete(limits = c("Pre-mature", "Mature"))+
   theme(axis.text.x = element_text(size = 14),
       axis.text.y = element_text(size = 14),
       axis.title =  element_text(size = 14),
       axis.title.y = element_text(size = 14),
       plot.title = element_text(size = 18))+
  scale_fill_manual(values = c("#BCBDDC","#807DBA" ,"#54278F"))+
  ylab("Number of TSS")
dev.off()
```

# k
```{r}
library(parallel)
library(stringr)
cell_de <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/celltype.theta_lineage.Rds") %>% 
  dplyr::mutate(Gene = unlist(str_split_fixed(.$TSS,"@",2))[,1],
                chr = unlist(str_split_fixed(unlist(str_split_fixed(.$TSS,"@",2))[,2],":",3))[,1],
                start = unlist(str_split_fixed(unlist(str_split_fixed(.$TSS,"@",2))[,2],":",3))[,2],
                strand = unlist(str_split_fixed(unlist(str_split_fixed(.$TSS,"@",2))[,2],":",3))[,3])


linea_psi_fil <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/cell.expre_tss.filter_lineage.Rds")

all <- lapply(split(linea_psi_fil,linea_psi_fil$Gene), function(df){
  tss_lst <- unique(df$TSS)
  types <- lapply(tss_lst, function(i){major(df,i)}) %>% do.call(c,.)
  
  return(types)
}) %>% do.call(c,.) %>% data.frame()
order <- c("Progenitor", "Mature")
color <- c('#F4BC58','#8DA0CB')
names(order) <- color

major <- function(Tss, tss){
  max(Tss$theta) -> H
  Tss[Tss$TSS== tss,]$theta -> t
  type <- ifelse(t == H, "Major", "Minor")
  type
}
linea_psi_fil[linea_psi_fil$Group=="Mature"]
lapply(split(linea_psi_fil, linea_psi_fil$Gene), function(df){
  M <- df[df$Group=="Mature",]
  P <- df[df$Group=="Progenitor",]
  if(dim(M)[1]==1 & dim(P)[1]>1 & length(intersect(P$TSS,M$TSS))>0){
    O_P <- cell_de[cell_de$Group=="Progenitor" & cell_de$Gene==unique(M$Gene),]
    tss_lst <- c(setdiff(P$TSS,M$TSS))
    tss_type <- lapply(tss_lst, function(i){major(O_P,i)}) %>% do.call(c,.)
    data.frame(type = rep("Lost",length(tss_type)),
               TSS = tss_lst,
               Gene = rep(unique(M$Gene),length(tss_type)),
               type2 = tss_type)
  }else if(dim(P)[1]==1 & dim(M)[1]>1& length(intersect(P$TSS,M$TSS))>0){
    O_M <- cell_de[cell_de$Group=="Mature" & cell_de$Gene==unique(M$Gene),]
    tss_lst <- c(setdiff(M$TSS,P$TSS))
    tss_type <- lapply(tss_lst, function(i){major(O_M,i)}) %>% do.call(c,.)
    data.frame(type = rep("Gained",length(tss_type)),
               TSS = tss_lst,
               Gene = rep(unique(M$Gene),length(tss_type)),
               type2 = tss_type)
  }else{
    return(NULL)
  }
}) %>% do.call(rbind,.) -> res
res[res$Gene=="Ctse",]
res %>% dplyr::group_by(type,type2) %>% dplyr::add_count(name = "num") %>% dplyr::select(type,type2,num) %>% .[!duplicated(.),] -> res
setDT(res)
res[, proportion := num/sum(num), by = type]


rbind.data.frame(res,
                 data.frame(type=c("All","All"),
                            type2=c("Major","Minor"),
                            num=c(1382,3227),
                            proportion=c(1382,3227)/sum(c(1382,3227)))) -> plot


```
```{r fig.height=4, fig.width=3}
plot$type <- factor(plot$type, levels = c("All","Gained","Lost","No-diff"))

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/0605_majorminor.pdf", height = 4, width = 3)


  ggplot(plot[plot$type!="No-diff",],
         aes(x=type,
             y=proportion,
             fill=type2,
             label=proportion))+
  geom_bar(stat = "identity",
           position ="fill",alpha=0.7,width = 0.8 )+
    theme_classic()+
    scale_fill_manual(values = c("#BEAED4" ,"#FDC086"))+
    ylab("Proportion of lost ATSs(%)")+
    scale_y_continuous()+#百分比加这个
    xlab(NULL)+
      theme(legend.position = "top",
          legend.title = element_blank(),
         legend.text = element_text(size = 12),
         axis.text = element_text(size = 12),
         axis.title  = element_text(size = 14)) +
    guides(fill = guide_legend( nrow = 1))
dev.off()

```


# l
```{r}

library(ChIPseeker)
library(org.Mm.eg.db)
linea_psi_fil <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/cell.expre_tss.filter_lineage.Rds")

txdb <- makeTxDbFromGFF("/mnt/raid66/Personal_data/xuzijie/ref/MusMus/release101/Mus_musculus.GRCm38.101.chr.gtf", format = "gtf")
obj_gr_lst <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/gain_lost.gr.Rds")


lapply(obj_gr_lst, function(peaks){
  peak_anno <- annotatePeak(
    peaks,
    addFlankGeneInfo = T,
    genomicAnnotationPriority = c("Promoter", "5UTR","Exon", "3UTR", "Intron",
    "Downstream", "Intergenic"),
    tssRegion = c(-50, 50),  # 定义启动子区域范围
    TxDb = txdb,  # 基因组注释数据库
    annoDb = "org.Mm.eg.db",      # 用于基因ID转换
    ignoreDownstream=F
    )
  peak_anno
}) -> peak_info


```


```{r}
col <- c("#a3c8dc","#92c875","#f19232","#a95728","#76509b")
names(col) <- c("Promoter","5' UTR","Exon","Intron","Other")
rbind.data.frame(data.frame(class=c(rep("All",5),rep("Gain",5),rep("Lost",5)),
                            type=c(rep(names(col), 3)),
                            per=c(75.1714678,8.2304527, 15.22634,0.6172839,0.7544582,
                                  56.0000000,7.3333333,35.33333,0.6666667,0.6666667,
                                  72.9508197,8.4699454,17.75956,0.2732240,0.5464481)), stringsAsFactors = F) -> plot_new
```

```{r fig.height=5, fig.width=4.5}
plot_new$type <- factor(plot_new$type, levels = names(col)) 
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/01DE/0607_regions.pdf", height = 5, width = 4.5)


  ggplot(plot_new,
         aes(x=class,
             y=per,
             fill=reorder(type,-per),
             label=per))+
  geom_bar(stat = "identity",
           position ="fill",alpha=0.9,width = 0.8 )+
    theme_classic()+
    scale_fill_manual(values = c(col))+
    ylab("Percentage of TSSs(%)")+
    scale_y_continuous()+#百分比加这个
    xlab(NULL)+
   theme(legend.position = "right",
          legend.title = element_blank(),
         legend.text = element_text(size = 12),
         axis.text = element_text(size = 12))+
    guides(fill = guide_legend( nrow = 6))
dev.off()

```