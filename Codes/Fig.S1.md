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
  library(dplyr)
  library(tidyr)
})
```
# a and b
```{r}
covid <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/covid19/V7/00quant/DF_rowRanges.Rds")
covid$type <- "COVID-19"
nsclc <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/00quant/NSCLC.DF_100.Rds")
nsclc@rowRanges$type <- "Lung cancer"
mbm <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds")
mbm@rowRanges$type <- "Mouse bone marrow"

rbind(data.frame(covid[,c("AUC","alpha1","type")]),
      data.frame(nsclc@rowRanges[,c("AUC","alpha1","type")]),
      data.frame(mbm@rowRanges[,c("AUC","alpha1","type")])) -> plot

setDT(plot)
plot[,.(AUC_mean=mean(AUC, na.rm=T), alpha_mean=mean(alpha1, na.rm=T),
        AUC_median=median(AUC, na.rm=T), alpha_median=median(alpha1, na.rm=T)), by=type]
```

```{r}
pair <- list(
             c("Mouse bone marrow","Lung cancer"),
             c("Lung cancer","COVID-19"),#**
             c("Mouse bone marrow","COVID-19")#*
             )#***

plot[,6:8] %>%reshape2::melt(.) -> plot2
mypal0 <- c("#e8d9de","#d49cad","#8f2041")
names(mypal0) <- c("Mouse bone marrow","Lung cancer","COVID-19")
plot2$type <- factor(plot2$type, levels = c("Mouse bone marrow","Lung cancer","COVID-19"))
```


```{r fig.height=4, fig.width=4}
ggplot(plot2[plot2$variable=="AUC",], aes(x=type, 
                 y=value,fill = type))+
  geom_violin(aes(fill = type),position=position_dodge(.9),color="white")+
  geom_boxplot(aes(fill = type),position=position_dodge(.9), width=0.2, outlier.shape = NA, color="black") +
  stat_boxplot(aes(fill = type),geom = "errorbar",#error bar
               position=position_dodge(.9),width=0.1,color="black")+
  stat_summary(fun = mean,#mean值连线
               geom = "point",
               aes(group = type, col = "red"))+
  scale_fill_manual(values=mypal0)+ 
  theme_bw()+ 
  theme(axis.ticks = element_blank(),
        legend.title = element_blank(),
        axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 12),
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12, face = "bold", angle = 0, vjust = 0.5),
        legend.position = "none")+
  labs(x = "", y = "β", title = "Sample degradation")+
  geom_signif(comparisons = pair,step_increase = 0.1,
        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
        test = t.test,y_position=c(1)) + ylim(0.4,1.3) -> p1

ggplot(plot2[plot2$variable=="alpha",], aes(x=type, 
                 y=value,fill = type))+
  geom_violin(aes(fill = type),position=position_dodge(.9),color="white")+
  geom_boxplot(aes(fill = type),position=position_dodge(.9), width=0.2, outlier.shape = NA, color="black") +
  stat_boxplot(aes(fill = type),geom = "errorbar",#error bar
               position=position_dodge(.9),width=0.1,color="black")+
  stat_summary(fun = mean,#mean值连线
               geom = "point",
               aes(group = type, col = "red"))+
  scale_fill_manual(values=mypal0)+ 
  theme_bw()+ 
  theme(axis.ticks = element_blank(),
        legend.title = element_blank(),
        axis.text.y=element_text(size = 10),
        axis.text.x=element_text(size = 12),
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12, face = "bold", angle = 0, vjust = 0.5),
        legend.position = "none")+
  labs(x = "", y = "α", title = "Sample degradation")+
  geom_signif(comparisons = pair,step_increase = 0.03,
        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
        test = t.test,y_position=c(1.3)) + ylim(0,2) -> p2
```


```{r fig.height=3.5, fig.width=8}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/alpha.pdf", height = 3.5, width = 8)
cowplot::plot_grid(p2,p1, nrow = 1)
dev.off()
```


# c
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

df <- as.data.frame(loci_mia0)
df2 <- df %>%
  select(gene_name, loci, alpha1)
gene_cmp <- df2  %>%
  group_by(gene_name)  %>%
  summarise(
    distal_alpha = alpha1[loci == 1][1],
    max_proximal_alpha = max(alpha1[loci > 1], na.rm = TRUE),
    proximal_higher = max_proximal_alpha > distal_alpha,
    n_loci = n(),
    .groups = "drop"
  )
n_total <- nrow(gene_cmp)
n_proximal_higher <- sum(gene_cmp$proximal_higher, na.rm = TRUE)
ratio_proximal_higher <- n_proximal_higher / n_total

n_proximal_higher
ratio_proximal_higher
```

```{r}


box_df <- gene_cmp[gene_cmp$proximal_higher=="TRUE" & gene_cmp$distal_alpha!=0,] |>
  select(gene_name = gene_name, distal_alpha, max_proximal_alpha) |>
  pivot_longer(
    cols = c(distal_alpha, max_proximal_alpha),
    names_to = "TSS_type",
    values_to = "alpha1"
  ) |>
  mutate(
    TSS_type = recode(TSS_type,
                      distal_alpha = "Distal ATS",
                      max_proximal_alpha = "Proximal ATS")
  )

pval <- wilcox.test(
  gene_cmp[gene_cmp$proximal_higher=="TRUE"& gene_cmp$distal_alpha!=0,,]$distal_alpha,
  gene_cmp[gene_cmp$proximal_higher=="TRUE"& gene_cmp$distal_alpha!=0,,]$max_proximal_alpha,
  paired = TRUE
)$p.value
```


```{r fig.height=4, fig.width=4}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/ats.alpha_proximaldistal4.pdf",height = 4,width = 4)

ggplot(box_df, aes(x = TSS_type, y = alpha1, fill = TSS_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha=0.7) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
  annotate("text", x = 1.5, y = max(box_df$alpha1, na.rm = TRUE),
           size=5,
           label = paste0("p = ", signif(pval, 3), ", N = 259")) +
  labs(
    x = NULL,
    y = "Degradation rate (alpha)",
    title = ""
  ) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()
```

# d
```{r}
gtfFile <- file.path("/mnt/raid61/Personal_data/tangchao/Document/gencode/mouse/release_M25/gencode.vM25.primary_assembly.annotation.sorted.gtf")
bams <- list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/bam/01bam","bam$",full.names = T)[2]
scats <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/00quant/TSSCDF_V8_withRD_0424.Rds")

gene_symbol <- "Mier1"
tss <- c(129887289,129887470)

scATS::Sashimi(object = scats, 
               bam = bams,
               xlimit = c(129886989,129887870),
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


# e
```{r}
c("#FC8D62","#66C2A5", "#ecd276","#da954f") -> col
names(col) <- c("scATS","SCAFE","TSSr","scTSS")

ref <- readRDS(paste0(outdir,"/quantify_benchmark/ref_hg38_tss.Rds"))

res <- readRDS(paste0(outdir,"/all_tool4_new.Rds"))

```

```{r fig.height=3, fig.width=4}
pdf(paste0(outdir, "/identificarion_4.pdf"), height = 3, width = 4)

ggplot(res, aes(x = class, y = value, fill = factor(type, levels = c("scATS","TSSr","SCAFE","scTSS")))) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7,alpha=0.9
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = col
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c("top")
  )+
  ylim(0,1)+
  guides(fill=guide_legend(nrow = 1))
dev.off()
```

# f
```{r}
lapply(Data, function(plot4){
  ct <- cor.test(plot4$psi, plot4$true_psi)
  
  lab <- paste(
    sprintf("R = %.2f", ct$estimate),
    paste0("N = ", length(unique(plot4$GENEID))),
    sep = "\n"
  )
  ggplot(plot4, aes(x = psi, y = true_psi)) +
      geom_hex(bins = 20) +
      geom_point(color = "#D48759", alpha = 0.15, size = 0.4) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      annotate(
          "text",
          x = -Inf, y = Inf,
          label = lab,
          hjust = -0.05,
          vjust = 1.1,
          size = 4
      ) +
      scale_fill_viridis_c() +
      theme_test() +
    theme(legend.position = "top")+
      labs(x = "scTSS", y = "Ground truth")->p
  p
  
}) -> p_lst


```
```{r fig.height=2.5, fig.width=9}
pdf(paste0(outdir, "ats_quantification.pdf"), height = 2.5, width = 9)
cowplot::plot_grid(p_lst, nrow = 1)
dev.off()
```


# g

```{r}
outdir <- "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/simulation/sequence_deepkinet/"
lapply(1:3, function(n){
  scats_lst <- readRDS( paste0(outdir,"/RD",n,"/all_20_scats_df_sc.Rds"))
  scats_lst$rep <- paste0("rep",n)
  scats_lst
}) %>% do.call(rbind,.) -> scats_lst
scats_lst$tool <- "scATS"
```


```{r}
load(paste0(outdir,"/SimObjects.RData"))

dir_name <- c("RD1","RD2", "RD3")
lapply(1:3, function(n){
  
  lapply(1:20, function(i){
    f <- paste0(outdir,dir_name[n],"/v",i,"/degradation_rate_single_cell_long.csv")
    print(file.exists(f))
    read.csv(f) %>% dplyr::mutate(version=i) ->df
    preset_list <- fromJSON(paste0(outdir,dir_name[n],"/v",i,"/genes_alpha.json"))
    preset_list2 <- do.call(
    rbind,
    lapply(names(preset_list), function(gene) {
      data.frame(
        gene  = gene,
        tss   = preset_list[[gene]]$tss_list,
        alpha = preset_list[[gene]]$alpha_vec,
        mean_alpha = mean(preset_list[[gene]]$alpha_vec),
        stringsAsFactors = FALSE
      )
    })
  )
    df$mean_alpha <- plyr::mapvalues(from = preset_list2$gene, to = preset_list2$mean_alpha, x = df$Gene, warn_missing = T) %>% as.numeric()
    df
  
}) %>% do.call(rbind,.) -> res_v1
  res_v1$rep <- paste0("rep",n)
  res_v1
}) %>% do.call(rbind,.) -> res_deepkinet
res_deepkinet$tool <- "DeepKINET"
saveRDS(res_deepkinet, file = paste0(outdir,"/res_deepkinet_sc.Rds"))


```

```{r}
res_deepkinet <- readRDS(paste0(outdir,"/res_deepkinet_sc.Rds"))
load(paste0(outdir,"/SimObjects.RData"))

res_deepkinet %>% dplyr::select(CellID,Gene,version,degradation_rate,mean_alpha,rep,tool) ->res_1
scats_lst %>% dplyr::select(alpha,gene,version,value,preset_alpha,rep,tool) ->res_2
colnames(res_1) <- c("cell","gene", "version", "alpha_quantify", "alpha_real", "rep","tool")
colnames(res_2) <- c("cell","gene", "version", "alpha_quantify", "alpha_real", "rep","tool")
rbind(res_1,res_2) -> res_plot

gene_with1tss <- names(which(table(gtf_tx0$gene_name)==1))
res_plot2 <- res_plot[res_plot$gene %in% gene_with1tss,]

res_plot2[!is.na(res_plot2$alpha_quantify)& res_plot2$alpha_quantify>0 &res_plot2$gene %in% gene_with1tss,] ->plot

bar_df1 <- plot %>%
  dplyr::group_by(cell,tool,gene,alpha_real) %>%
  dplyr::mutate(
    mean_alpha = mean(alpha_quantify,na.rm=T)
  )%>% dplyr::ungroup() %>%
  dplyr::select(alpha_real,tool,mean_alpha,gene) %>%
  dplyr::distinct()
saveRDS(bar_df1,paste0(outdir,"/plot/RD/sample_sc.Rds"))

```


```{r fig.height=2.5, fig.width=8}
pdf(paste0(outdir,"/plot/RD/sample_sc.pdf"), height = 2.5, width = 8)

ggplot(bar_df1, aes(
  x = factor(alpha_real), 
  y = mean_alpha,
  fill = factor(tool, levels = rev(c("DeepKINET","scATS")))
)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.color = NA,
    alpha = 0.6
  ) +
  stat_boxplot(geom = "errorbar",#error bar
               position=position_dodge(width=0.75))+
  labs(
    x = "Ground truth ",
    y = "Estimated degradation rate"
  ) + ylim(0,3)+
  theme_test()+
  theme(legend.position = c(0.15, 0.92),
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  scale_color_manual(values = color)+
  scale_fill_manual(values = color)+
  guides(color=guide_legend(nrow=1),
         fill=guide_legend(nrow=1))
dev.off()

```


# h
```{r}
outdir <- "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/simulation/sequence_deepkinet/"
lapply(1:3, function(n){
  scats_lst <- readRDS( paste0(outdir,"/RD",n,"/all_20_scats_df_sc.Rds"))
  scats_lst$rep <- paste0("rep",n)
  scats_lst
}) %>% do.call(rbind,.) -> scats_lst
scats_lst$tool <- "scATS"
```

```{r}
res_deepkinet <- readRDS(paste0(outdir,"/res_deepkinet.Rds"))
load(paste0(outdir,"/SimObjects.RData"))
load(paste0(outdir,"/gene_with_1or2_tss.RData"))

res_deepkinet %>% dplyr::select(gene,version,degradation_rate_mean,mean_alpha,rep,tool) ->res_1
scats_lst %>% dplyr::select(gene_name,version,alpha,preset_alpha,rep,tool) ->res_2
colnames(res_1) <- c("gene", "version", "alpha_quantify", "alpha_real", "rep","tool")
colnames(res_2) <- c("gene", "version", "alpha_quantify", "alpha_real", "rep","tool")
rbind(res_1,res_2) -> res_plot

gene_with1tss <- names(which(table(gtf_tx0$gene_name)==1))
res_plot2 <- res_plot[res_plot$gene %in% gene_with1tss,]
cor_df0 <- res_plot2 %>%
  group_by(tool,rep) %>%
  summarise(
    cor_alpha = cor(alpha_quantify, alpha_real)
  )%>% ungroup()

cor_df1 <- res_plot2 %>%
  group_by(tool,rep,gene) %>%
  summarise(
    cor_alpha = cor(alpha_quantify, alpha_real)
  )%>% ungroup()

```

```{r fig.height=3, fig.width=4}
color <- rev(c("#63a4c9","#FC8D62"))
names(color) <- rev(c("DeepKINET","scATS"))
cor_df1[!is.na(cor_df1$cor_alpha)& cor_df1$cor_alpha>0,] %>%
  dplyr::group_by(tool) %>%
  dplyr::summarise(mean_v=mean(cor_alpha, na.rm=T))

pdf(paste0(outdir,"/plot/RD/rep3.pdf"), height = 3, width = 4)
ggplot(cor_df1[!is.na(cor_df1$cor_alpha)& cor_df1$cor_alpha>0,], 
       aes(x=factor(tool, levels = rev(c("DeepKINET","scATS"))),
           y=cor_alpha, color=rep))+
  geom_violin(alpha=0.8)+
  # geom_boxplot(alpha=0.8)+
  theme_bw()+
  labs(x="",y="Correlation")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  scale_color_brewer(palette = "Dark2")+
  guides(color=guide_legend(nrow=1))
dev.off()
```