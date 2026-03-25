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

# e
```{bash}

python /mnt/raid61/Personal_data/xuzijie/task/script/geneBody_coverage.py \
-i /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/bam/01bam/merge.bam \
-r /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/gencode.vM25.primary_assembly.annotation.exon.bed12 \
-o total

python /mnt/raid61/Personal_data/xuzijie/task/script/geneBody_coverage.py \
-i /mnt/raid64/Covid19_Gravida/cellranger/CG/outs/possorted_genome_bam.bam \
-r /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/covid19/Homo_sapiens.GRCh38.101.exon.bed12 \
-o covid19

python /mnt/raid61/Personal_data/xuzijie/task/script/geneBody_coverage.py \
-i /mnt/raid61/Personaln_data/xuzijie/task/07ATS/00data/thesis/NSCLC/bam//NSCLC.bam \
-r /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/covid19/Homo_sapiens.GRCh38.101.exon.bed12 \
-o nsclc

python /mnt/raid61/Personal_data/xuzijie/task/script/geneBody_coverage.py \
-i /mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/Smart-seq/bams/Smartseq_BM42.bam \
-r /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/Mus_musculus.GRCm38.101.exon.bed12 \
-o Smartseq

python /mnt/raid61/Personal_data/xuzijie/task/script/geneBody_coverage.py \
-i /mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/3scrnaseq/bams/sc_NGScDNA_C57_F_2m_c-kit+.bam \
-r /mnt/raid61/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/Mus_musculus.GRCm38.101.exon.bed12 \
-o sc_NGScDNA_C57_F_2m_c-kit+
  
```

```{r}

plot_coverage <- function(files,name,color){
  lapply(files,function(x){
    tab <- fread(x,sep = '\t', header = F)
    tab <- tab[2,]
  }) %>% do.call(rbind,.) -> geneBodyCoverage
  
  geneBodyCoverage2<- apply(geneBodyCoverage[,-1],1,function(x){
    (x - min(x)) / (max(x) - min(x))   
  })
  colnames(geneBodyCoverage2) <- name
  geneBodyCoverage2 <- cbind(geneBodyCoverage2,percent=1:100)
  geneBodyCoverage2 <- melt(geneBodyCoverage2,id.var="percent")
  geneBodyCoverage2$Var1<-rep(1:100,length(unique(geneBodyCoverage2$Var2)))
  geneBodyCoverage2 <-
    geneBodyCoverage2[geneBodyCoverage2$Var2 != 'percent', ]
  
  
  unique(geneBodyCoverage2$Var2)-> names(color)
  p <-
    ggplot(geneBodyCoverage2,
           aes(
             x = Var1,
             y = value,
             group = Var2,
             colour = Var2
           )) + geom_line() +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(),
      title = element_text(
        size = 15,
        hjust = .5,
        vjust = .5
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      plot.title = element_text(hjust = 0.5,size = 18),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(x = "5' -> 3'", y = "Coverage")+
    scale_color_manual(values = c(color))
  
  return(p)
}
```

```{r}
c(
  "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/00sequence/total.geneBodyCoverage.txt",
  "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/00sequence/nsclc.geneBodyCoverage.txt",
  "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/00sequence/covid19.geneBodyCoverage.txt",
  "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V6/00quant/coverage/smart/Smartseq.geneBodyCoverage.txt",
  "/mnt/raid66/Personal_data/zhangdan/Isoform_atlas/01_data_ingest/11_RSEQC/gene_body/sc_NGScDNA_C57_F_2m_c-kit+.geneBodyCoverage.txt" ) -> file

mypal<-c("#e2bdc8",
         "#d49cad",
         "#8f2041",
         "#a7a88e",
         "#255144")
p <- plot_coverage(file, c("5scRNA-seq(mm_BM)",
                           "5scRNA-seq(NSCLC)",
                           "5scRNA-seq(Covid)",
                           "Smart-seq(mm_BM)",
                           "3scRNA_seq(mm_BM)"),
                   mypal)

ggsave(p, file = paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V6/00quant/coverage/genebody.pdf"), width = 6, height = 3)
```


```{r}
files <- file
name <- c("5scRNA-seq(mm_BM)","5scRNA-seq(NSCLC)","5scRNA-seq(Covid)",
                           "Smart-seq(mm_BM)","3scRNA_seq(mm_BM)")
color <- c("#1587AF", "#EAA61F","#93CC82","#7D4FC6","#B1B1B2")
lapply(files,function(x){
  tab <- fread(x,sep = '\t', header = F)
  tab <- tab[2,]
}) %>% do.call(rbind,.) -> geneBodyCoverage

geneBodyCoverage2<- apply(geneBodyCoverage[,-1],1,function(x){
})
colnames(geneBodyCoverage2) <- name
geneBodyCoverage2 <- cbind(geneBodyCoverage2,percent=1:100)
geneBodyCoverage2 <- melt(geneBodyCoverage2,id.var="percent")
geneBodyCoverage2$Var1<-rep(1:100,length(unique(geneBodyCoverage2$Var2)))
geneBodyCoverage2 <-
  geneBodyCoverage2[geneBodyCoverage2$Var2 != 'percent', ]


unique(geneBodyCoverage2$Var2)-> names(color)
```


```{r fig.width = 6, fig.height = 2}
ggplot(geneBodyCoverage2[geneBodyCoverage2$Var1%in% 1:10,],
       aes(
         x = Var1,
         y = value,
         group = Var2,
         colour = Var2
       )) + geom_line() +
theme(
  panel.background = element_blank(),
  axis.line = element_line(),
  title = element_text(
    size = 15,
    hjust = .5,
    vjust = .5
  ),
  legend.title = element_blank(),
  legend.text = element_text(size = 13),
  axis.text = element_text(size = 13),
  axis.title = element_text(size = 15),
  plot.title = element_text(hjust = 0.5,size = 18),
  plot.subtitle = element_text(hjust = 0.5),
  legend.position = "none"
) +
labs(x = "5' -> 3'", y = "Coverage")+
scale_color_manual(values = c(color))
```


```{r fig.width = 6, fig.height = 3}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/coverage.genebody_firstexon.pdf", height = 3, width = 6)
p
dev.off()
```


# f
```{r}
gtfFile <- file.path("/mnt/raid61/Personal_data/tangchao/Document/gencode/mouse/release_M25/gencode.vM25.primary_assembly.annotation.sorted.gtf")
bams <- list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/bam/01bam","bam$",full.names = T)[2]
scats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds")

gene_symbol <- "Chchd2"
tss <- c(129887289,129887470)

scATS::Sashimi(object = scats, 
               bam = bams,
               xlimit = c(129886989:129887870),
               gtf = gtfFile, 
               gene = gene_symbol, 
               TSS = tss, #等于line
               free_y = T,#是否scale
               base_size = 12, #read部分字体大小
               rel_height=0.5, #注释/read ，小于1 read部分比例更大
               line.color = "black",
               line.type = 3) -> p 
p
```

# g
```{r}
scats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds")
scats@rowRanges %>%data.frame()->all_loci

all_loci2 <- all_loci[!is.na(all_loci$alpha1),]
lst=seq(0, ceiling(max(all_loci2$alpha1)) , 0.2)

lapply(1:(length((lst))-1) ,function(x){
  df0 <- all_loci2[all_loci2$alpha1 < lst[x+1] & all_loci2$alpha1 >= lst[x],]
  df_precision <- data.frame(value=dim(df0)[1] * 100/dim(all_loci2)[1],
                             num = dim(df0)[1],
                             cutoff=paste0("[",lst[x],",",lst[x+1],")"),
                             class="All") 
  df_precision
}) %>% do.call(rbind,.) -> plot_all
plot_all$label <- paste0("(",plot_all$num,") ", plot_all$cutoff)
plot_all$cutoff2 <- ifelse(plot_all$cutoff %in% c(as.character(plot_all$cutoff)[1:10]),as.character(plot_all$cutoff),"[2,+]" )

```


```{r fig.width=4,fig.height=3}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/ats.cutoff_alpha.pdf",height = 3,width = 4)
ggplot(plot_all, aes(x = cutoff2, y = value, color = class,group=class))+
    geom_point(size = 2, alpha = 0.9) +
    geom_line(size = 1, alpha = 0.9) +
    labs(x = "Cutoff of alpha", y = "Proportion of all ATSs(%)") +
    theme_light()+
  scale_color_manual(values = "#3C5488E5")+
    theme(#legend.position=c(.2,0.95), 
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      axis.title.x = element_text(hjust = 0.5,size = 14),
      axis.title.y =  element_text(hjust = 0.5,size = 14),
      axis.text.y= element_text(color="black", size=12),
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        legend.position = "none"
    ) 
dev.off()
```

# h and i

```{r}
scats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds")
scats@rowRanges %>%data.frame()->all_loci


mclapply(split(all_loci,all_loci$gene_name), function(df){
  if(length(unique(df$start)) > 1){ # 去掉只有一个ats的基因
  
    if(df$strand == "+"){
      df <- df[order(df$start,decreasing = F),]
      df$loci <- c(1:dim(df)[1])
      return(df)
    }else{
      df <- df[order(df$start,decreasing = T),]
      df$loci <- c(1:dim(df)[1])
      return(df)
    }
    
  }
},mc.cores = 20) %>% do.call(rbind,.) -> loci_mia0
loci_mia0$loci <- ifelse(loci_mia0$loci>=4, ">=4",loci_mia0$loci)

loci_mia0 %>% dplyr::group_by(loci) %>% dplyr::add_count(name = "n_loci")-> loci_mia0
loci_mia0$loci2 <- paste0(loci_mia0$loci,"(", loci_mia0$n_loci,")")

```


```{r fig.height=7, fig.width=4}
loci_mia0[,c("alpha1","AUC","loci2")] %>% reshape2::melt() -> plot_alpha

plot_alpha[plot_alpha$variable=="AUC",] -> plot2
ggplot(plot2, aes(x = factor(loci2, levels = c(lo)), y = value,fill=variable)) +
  geom_violin(position = position_dodge(width = 0.9))+
  geom_boxplot( outlier.shape = NA,position = position_dodge(width = 0.5), width=0.3) +
  stat_summary(fun = mean,#mean值连线
               geom = "point",
               aes(group = variable),
               color="red") +
stat_summary(fun = mean,#mean值连线
               geom = "line",
               aes(group = variable)) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c("#e6eef3")) +
  labs(y = "β", x = "ATS locations (distal to proximal)")+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
        strip.background =  element_rect(colour="black",
                                      fill="white"),
        strip.text.y = element_text(size = 10)) +
  geom_signif(comparisons = group,step_increase = 0.1,
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    test = wilcox.test) -> p1

plot_alpha[plot_alpha$variable=="alpha1",] -> plot2
ggplot(plot2, aes(x = factor(loci2, levels = c(lo)), y = value,fill=variable)) +
  geom_violin(position = position_dodge(width = 0.9))+
  geom_boxplot(outlier.shape = NA,position = position_dodge(width = 0.5), width=0.3) +
  stat_summary(fun = mean,#mean值连线
               geom = "point",
               aes(group = variable),
               color="red") +
stat_summary(fun = mean,#mean值连线
               geom = "line",
               aes(group = variable)) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c("#84a8c2")) +
  labs(y = "α", x = "ATSs per gene")+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background =  element_rect(colour="black",
                                      fill="white"),
        strip.text.y = element_text(size = 10)) +
  geom_signif(comparisons = group,step_increase = 0.05,
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    na.rm=T,y_position = 1.5,
    test = wilcox.test) + ylim(0,2)-> p2
```


```{r fig.height=7, fig.width=4}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/ats.alpha_proximaldistal3.pdf",height = 7,width = 4)
cowplot::plot_grid(p2,p1,nrow = 2)
dev.off()
```

# j
```{r fig.height=3, fig.width=3}

scats <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds") #
data.frame(scats@rowRanges) -> locis
locis <- locis[!is.na(locis$alpha1),]

locis2 <- locis[locis$gene_name %in% names(which(table(locis$gene_name)>1)),]
locis2[, c("PSI", "theta", "alpha1")] -> plot
plot[plot$alpha1>0.8,] -> plot2
plot <- plot[order(plot$alpha1, decreasing = F),]

plot$residual <- plot$theta - plot$PSI

sd_2 <- 2 * sd(plot$residual, na.rm = TRUE)
plot$outlier <- ifelse(abs(plot$residual) > sd_2, "Outlier", "Normal")
outlier_count <- sum(abs(plot$residual) > sd_2, na.rm = TRUE)
sd_po <- dim(plot[plot$residual>0 & plot$outlier=="Outlier",])[1]
sd_neg <- dim(plot[plot$residual<0 & plot$outlier=="Outlier",])[1]
```


```{r fig.height=3, fig.width=4}
locis$labels <- 0
locis[locis$gene_name %in% c("Chchd2","Hsd17b10"),"labels"] <- 1
locis[, c("PSI", "theta", "alpha1","labels","TSS")] -> plot

plot <- plot[order(plot$alpha1, decreasing = F),]
plot$residual <- plot$theta - plot$PSI
sd_2 <- 2 * sd(plot$residual, na.rm = TRUE)
plot$outlier <- ifelse(abs(plot$residual) > sd_2, "Outlier", "Normal")
outlier_count <- sum(abs(plot$residual) > sd_2, na.rm = TRUE)
sd_po <- dim(plot[plot$residual>0 & plot$outlier=="Outlier",])[1]
sd_neg <- dim(plot[plot$residual<0 & plot$outlier=="Outlier",])[1]

subset(plot,labels==1)
```


```{r fig.height=3, fig.width=4.5}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/0528.cor.psi.theta_all.pdf", width = 4.5, height = 3)
showtext::showtext.begin()
ggplot(plot, aes(x = PSI, y = theta, color = alpha1)) +

  # 点图层（按alpha1分层）
  geom_point(
    data = subset(plot, alpha1 <= 0.8),
    size = 1, 
    alpha = 0.6  # 增加透明度
  ) +
  geom_point(
    data = subset(plot, alpha1 > 0.8),
    size = 2,
    shape = 21,  # 添加边框
    fill = "red", 
    alpha=0.8,
    stroke = 0.5  # 边框粗细
  ) +
  
  # 颜色标度（重点修改部分）
  scale_color_steps2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(0, 1),  # 强制颜色范围0-1
    breaks = seq(0, 1, 0.2),  # 明确刻度间隔
    guide = guide_colorbar(  # 优化图例
      barwidth = unit(.5, "cm"),
      title.position = "top"
    )
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "gray40", linewidth = 0.5) +
  geom_abline(intercept = sd_2, slope = 1, linetype = "dotted", 
              color = "gray20", linewidth = 0.5) +
  geom_abline(intercept = -sd_2, slope = 1, linetype = "dotted", 
              color = "gray20", linewidth = 0.5) +
  annotate("text", x = 0.95, y = 0.9 , 
           label = "y=x", color = "gray20", size = 4) +
  # 标注2SD范围
  annotate("text", x = max(plot$PSI)*0.9-0.1, y = max(plot$PSI)*0.9 + sd_2+0.07 , 
           label = paste0("+2SD"," (N=",sd_po,")"), color = "gray20", size = 3) +
  annotate("text", x = max(plot$PSI)*0.9, y = max(plot$PSI)*0.9 - sd_2 -0.3, 
           label = paste0("-2SD"," (N=",sd_neg,")"), color = "gray20", size = 3)  +
  geom_text(data=subset(plot,labels==1), aes(x =PSI ,y =theta ,label=TSS),size=3,color="black",
            direction="y",hjust=0) +
  geom_point(
    data = subset(plot,labels==1),
    size = 2,
    shape = 21,  # 添加边框
    fill = "red", 
    color="black",
    alpha=0.8,
    stroke = 0.5  # 边框粗细
  )+
  geom_vline(xintercept = 0.9)+
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., paste0("`N = ", length(plot$PSI), "`"), sep = "~`,`~")),
    label.x = 0, 
    label.y = 1.2, 
    color = 'black', 
    size = 4
  ) +
  # 主题与标签
  labs(
    x = "Φ", 
    y = "θ", 
    title = "Genes with all TSSs",
    color = "α")+  # 数学符号标注
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14,face = "bold"),
    axis.title.y = element_text(size = 14,face = "bold",angle = 0, vjust = 0.5),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )
showtext::showtext.end()
dev.off()
```


# k
```{r}
scats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds")
scats@rowRanges %>%data.frame()->all_loci


mclapply(split(all_loci,all_loci$gene_name), function(df){
  if(length(unique(df$start)) > 1){ # 去掉只有一个ats的基因
  
    if(df$strand == "+"){
      df <- df[order(df$start,decreasing = F),]
      df$loci <- c(1:dim(df)[1])
      return(df)
    }else{
      df <- df[order(df$start,decreasing = T),]
      df$loci <- c(1:dim(df)[1])
      return(df)
    }
    
  }
},mc.cores = 20) %>% do.call(rbind,.) -> loci_mia0
loci_mia0$loci <- ifelse(loci_mia0$loci>=4, ">=4",loci_mia0$loci)

loci_mia0 %>% dplyr::group_by(loci) %>% dplyr::add_count(name = "n_loci")-> loci_mia0
loci_mia0$loci2 <- paste0(loci_mia0$loci,"(", loci_mia0$n_loci,")")

```

```{r}
library(ggsignif)
unique(loci_mia0$loci2) -> lo
group <- list(
              c(lo[1], lo[2]),
              c(lo[1], lo[3]),
              c(lo[1], lo[4]))
loci_mia0$Apices_del <- loci_mia0$ApicesY - loci_mia0$ApicesX
loci_mia0[,c("PSI","theta","loci2")] %>% reshape2::melt() -> plot_alpha
```
```{r}
loci_mia0[loci_mia0$gene_id %in% names(which(table(loci_mia0$gene_id)==2)),] -> loci_mia2

```


```{r fig.width=4,fig.height=3}
levels(plot_alpha$variable) <- c("theta","PSI")
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/0414_distal_psi_theta2.pdf", height = 3, width = 4)

ggplot(plot_alpha, aes(x = factor(loci2, levels = c(lo)), y = value)) +
  geom_violin(alpha=0.9,position = position_dodge(width = 0.7),aes(fill=variable))+
  geom_boxplot(alpha=0.9,aes(fill=variable), outlier.shape = NA,position = position_dodge(width = 0.7), width=0.2) +
  stat_summary(fun = mean,#mean值连线
               geom = "point",
               aes(group = variable,fill=variable),position = position_dodge(width = 0.7),
               color="red") +
stat_summary(fun = mean,#mean值连线
               geom = "line",position = position_dodge(width = 0.7),
               aes(group = variable,color=variable)) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c("#399765", "#b3dba8")) +
  scale_color_manual(values = c("#399765", "#b3dba8")) +
  labs(y = "PSI", x = "ATS location (distal to proximal)") +
  theme( legend.title = element_blank(),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12),
        legend.position = c(0.8, 0.8),
        strip.background =  element_rect(colour="black",
                                      fill="white"),
        strip.text.y = element_text(size = 10)) +
  stat_compare_means(aes(group=variable), label = "p.signif")
dev.off()
```
