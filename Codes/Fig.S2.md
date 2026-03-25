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

# b

```{r}
tss_df2 <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATS_obj_5_0105.Rds")

gtf_dis_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/distance/all5_gtf.Rds")
reftss_dis_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/distance/all5_reftss.Rds")

limit <- 100
levels <- c("scATS","SCAFE","CamoTSS", "TSSr", "scTSS")
lapply(1:5, function(i){
  
  tc_gtf_df <- gtf_dis_lst[[i]][,c("site","strand","distance")]
  tc_gtf_df$distance[tc_gtf_df$distance >= limit] <- limit
  tc_gtf_df$distance[tc_gtf_df$distance <= -limit] <- -limit
  tc_gtf_df$distance <- abs(tc_gtf_df$distance)
  tc_gtf_df$class <- levels[i]
  tc_gtf_df
}) -> gtf_lst

limit <- 50
lapply(1:5, function(i){

  tc_reftss_df <- reftss_dis_lst[[i]][,c("site","strand","distance")]
  tc_reftss_df$distance[tc_reftss_df$distance >= limit] <- limit
  tc_reftss_df$distance[tc_reftss_df$distance <= -limit] <- -limit
  tc_reftss_df$distance <- abs(tc_reftss_df$distance)
  tc_reftss_df$class <- levels[i]
  tc_reftss_df
}) -> reftss_lst
```

```{r fig.height=4, fig.width=4}

my1 <- list(c("scATS","scTSS"),
           c("scATS","CamoTSS"),
           c("scATS","SCAFE"),
           c("scATS","TSSr"))
my2 <- list(c("scATS","SCAFE"),
           c("scATS","CamoTSS"),
           c("scATS","scTSS"),
           c("scATS","TSSr"))
library(ggpubr)
library(gridExtra)
plot_fig1 <- function(plot, title, my){
  ggplot(plot,aes(x=distance,color=class))+
    geom_density()+
    scale_color_manual(values = col, name='')+
    theme_classic()+
      theme(legend.position = "none",
            axis.text.y = element_text(size = 12),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
         axis.title = element_text(size = 14))+ylab("Density")+
      xlab('') ->p1
  d <- max(plot$distance)
  ggplot(plot,aes(x=class,y=distance,color=class))+
    geom_boxplot(position=position_dodge(.9), width=0.5,
                 outlier.shape = NA)+
      stat_boxplot(geom = "errorbar",#error bar
                   position=position_dodge(.9),width=0.5)+
    scale_color_manual(values = col, name='')+
    theme_classic()+
      geom_signif(comparisons = my,step_increase = 0.1,
        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
        test = wilcox.test,y_position=c(d,d+1,d+1.5, d+2))+
      theme(legend.position = "none",
            axis.text = element_text(size = 12),
         axis.title = element_text(size = 14))+xlab("")+
      ylab(title)+coord_flip() -> p2
  
  combined_plot <- p1 / p2 + patchwork::plot_layout(heights = c(3, 1.5))
  combined_plot
} 
```


```{r fig.height=3, fig.width=3}
or1 <- c("scATS","scTSS","CamoTSS","SCAFE","TSSr")
or2 <- c("scATS","SCAFE","CamoTSS","scTSS","TSSr")

pdf('/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/distance/distance_densities_1206.pdf', height = 3, width = 3)
rbind(gtf_lst[[1]], gtf_lst[[5]], gtf_lst[[3]], gtf_lst[[2]], gtf_lst[[4]])-> plot
plot$class <- factor(plot$class, levels = rev(or1))
plot_fig1(plot,title="Distance to GENCODE TSSs",my1)

rbind(reftss_lst[[1]], reftss_lst[[2]], reftss_lst[[3]], reftss_lst[[4]], reftss_lst[[5]])-> plot
plot$class <- factor(plot$class, levels = rev(or2))
plot_fig1(plot,title="Distance to refTSS  TSSs",my2)

dev.off()
```


# c
```{r}
ALL_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATSRD_SCAFE_camotss_obj8_0106.Rds")

lapply(seq_along(ALL_lst),function(i){
  data.frame(ALL_lst[[i]])->df
  df$peak <- paste0("peak_",c(1:dim(df)[1]))
  
  df$end <- ifelse(df$strand=="+",df$end +1,df$end)
  df$start <- ifelse(df$strand=="-",df$start -1,df$start)
  
  write.table(df[,c("peak","seqnames","start","end","strand")],file = paste0("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V7/02accuracy/homer_1123/",names(ALL_lst)[i],".bed"),col.names = F,row.names = F,quote = F,sep = "\t")
})
```
```{bash}
for x in `cat data.txt`
do
echo $x
for i in `cat Initiator.table.txt`
do
echo $i
/mnt/data3/xuzijie/myconda3/envs/python3.8/share/homer/bin/findMotifsGenome.pl $x.bed mm10 ./ -find /mnt/data3/xuzijie/myconda3/envs/python3.8/share/homer/data/knownTFs/motifs/$i -size -50,50 > $x/$i.txt
done
done


```

```{r}
library(dplyr)
bed <- lapply(list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V7/02accuracy/homer_1123",".bed",full.names = T)[1:8],function(i){
  print(i)
  read.delim2(i,header = F,stringsAsFactors = F)->df
  colnames(df) <- c("peak","seqnames","start","end","strand")
  df$end <- as.numeric(df$end)
  df$start <- as.numeric(df$start)
  df$end <- ifelse(df$strand=="+",df$end -1,df$end)
  df$start <- ifelse(df$strand=="-",df$start +1,df$start)
  df$site <- paste0(df$seqnames,":",df$start,":",df$strand)
  df
})

homer <- lapply(list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V7/02accuracy/homer_1123",full.names = T),function(dir){
  do.call(rbind,lapply(list.files(dir,"motif.txt$",full.names = T),function(i){read.delim2(i,header = T,stringsAsFactors = F)}))->df
  df
})
names(homer) <-  basename(list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V7/02accuracy/homer_1123",full.names = T))


name <- names(homer) 
lapply(name,function(i){
  print(i)
  df <- merge(bed[[i]],homer[[i]],by.x = "peak",by.y = "PositionID")
  if(i=="CAGE"){
    df1 <- df[df$site %in% ALL_lst[[i]]$site,]
  }else{
    df1 <- df[df$site %in% as.character(ALL_lst[[i]]),]
  }
  df1 %>% dplyr::select(peak,site,Offset,Motif.Name,MotifScore)->df2
  df2$class <- i
  df2
})->peak_data
names(peak_data) <- name
peak_data$CAGE <- peak_data$CAGE[peak_data$CAGE$site %in% ALL_lst$CAGE$site,]


lapply(names(peak_data),function(i){
  data.frame(type = i,
             prop = round(length(unique(peak_data[[i]]$peak)) / length(unique(bed[[i]]$peak))* 100),
             num = length(unique(peak_data[[i]]$peak)))
  
})%>% do.call(rbind,.)->data_coun
data_coun

```


```{r fig.height=3.2,fig.width=3}
library(ggplot2)
library(ggpubr)

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/homer3/all8_3.pdf", height = 3.2, width = 3)
ggplot(data_coun,
       aes(x=reorder(type, -prop),
           y=prop,
           fill=type))+
geom_bar(stat = "identity",
         position ="dodge",width=0.7,alpha=0.7)+
  # theme_classic()+
  geom_text(aes(label=num),
            position=position_dodge(.9),
            size=4,hjust=.5, vjust=1)+
  scale_fill_manual(values = col, name='')+
  theme_bw()+
  labs(y="Percentage of TSSs (%)", x="", title = "Promoter motifs (±50bp)")+
 theme(legend.position = "none",
        axis.title = element_text(size = 14),
       axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
       axis.text.y = element_text(size = 12))  
dev.off()
```

# d
```{r}

library(magrittr)
library(plyr)
library(stringr)
library(GenomicRanges)

ALL_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATSRD_SCAFE_camotss_obj8_0106.Rds")

CpG_g <- readRDS("/mnt/raid66/Personal_data/xuzijie/ref/MusMus/UCSC/mm10.CpG.Rds")

file <- "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/CpG/1123/"
lapply(c(1:7), function(i){
  get_dis_pices_core(ALL_lst[[i]], CpG_g, path= paste0(file,names(ALL_lst)[i],"/"))
})


dir <- list.files(file,full.names = T)
lapply(dir,function(i){
  lapply(list.files(i,pattern = "Rds",full.names = T),function(y){
    readRDS(y)
  })%>% do.call(rbind,.)->df
  df
})->data_lst
names(data_lst) <- basename(dir)
data_lst2 <- data_lst
```


```{r}
library(dplyr)
library(data.table)
seq_limit <- 2500

data_plot2 <- lapply(c(1:8),function(i){

  df <- data_lst2[[i]]
  df$distance <- as.numeric(df$distance)
  df2 <- df[df$distance >= -2000 & df$distance < 500,]
  data.frame(type =  names(data_lst2)[i],
             num = dim(df2)[1],
             prop = round(dim(df2)[1]/dim(df)[1] * 100),
             count_mean = mean(na.omit(as.numeric(df2$count))),
             count_median = median(na.omit(as.numeric(df2$count))))

})%>% do.call(rbind,.)
data_plot2

```


```{r fig.height=3.2,fig.width=3}

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/CpG/all8_3.pdf", height = 3.2, width = 3)
ggplot(data_plot2,
       aes(x=reorder(type, -prop),
           y=prop,
           fill=type))+
geom_bar(stat = "identity",
         position ="dodge",width=0.7,alpha=0.7)+
  theme_bw()+
  geom_text(aes(label=num),
            position=position_dodge(.9),
            size=4,hjust=.5, vjust=1)+
  scale_fill_manual(values = col, name='')+
  labs(y="Percentage of TSSs (%)", x="", title = "CpG island (-2kb ~ +0.5kb)")+
 theme(legend.position = "none",
        axis.title = element_text(size = 14),
       axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
       axis.text.y = element_text(size = 12))+
  ylim(0,100)
dev.off()
```

# e
```{r}
ALL_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATSRD_SCAFE_camotss_obj8_0106.Rds")

bw <- list.files("/mnt/raid64/ATS/rawdata/ENCODE/mm10/MPP","bigWig$",full.names = T,recursive = T)

lapply(bw,function(b){
  print(b)
  lapply(c(1:8),function(i){
    print(i)
    res = ALL_lst[[i]]
    lapply(c(1:2),function(n){
      phast <- get_hist(g_rds = split(res, strand(res))[[n]],
                        bw_path = b,
                        seq = 2500,
                        label = paste0(names(ALL_lst)[i],"_",names(split(res, strand(res)))[n])) %>%
        mutate(Histone = gsub(".bigWig","",str_split_fixed(b,"/",12)[,9]))
        
      phast
})%>% do.call(rbind,.)->a
  a
  })%>% do.call(rbind,.)->b
  b
})%>% do.call(rbind,.)->tc_lst

tc_lst$type1 <- mapply(function(x) x[1], strsplit(tc_lst$ss, "_"))
tc_lst$type2 <- mapply(function(x) x[2], strsplit(tc_lst$ss, "_"))
tc_lst[tc_lst$type2=="-",]$Pos <- -tc_lst[tc_lst$type2=="-",]$Pos
plot <- tc_lst[,.(mean = mean(y)), by=.(Histone,type1,Pos)]

```

```{r fig.width = 7,fig.height = 3}
level <- c("scATS","SCAFE","CamoTSS","TSSr","scTSS","GENCODE","CAGE","Random")
plot$Histone <- factor(plot$Histone, levels = c( "H3K4me3","H3K27ac","H3K4me1"))
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/histone/1125//all.8.pdf",height = 3,width = 7)

ggplot(plot, aes(x = Pos, y = mean,colour = factor(type1, levels = level))) + 
    geom_line(size = 0.7,alpha=.7) + facet_wrap(.~Histone,scales = "free_y")+
  scale_color_manual(values =col)+
  ylab("Average coverage")+xlab("Distance to TSS")+
    scale_x_continuous(breaks = seq(-2500,2500,1250))+
    guides(shape = FALSE, colour = guide_legend(override.aes = list(size = 3), nrow = 1)) + 
    theme_classic()+
    theme(axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          axis.title =  element_text(size = 14),
          legend.text=element_text(size=14),
          plot.title = element_text (hjust = 0.5),
          legend.title = element_blank(),
          strip.text = element_text(size = 14),
          legend.position = "top")
dev.off()
```
# f
```{r}

ALL_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATSRD_SCAFE_camotss_obj8_0106.Rds")

bw <- "/mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/ENCODE/mouse/ATAC/HSC/ENCFF351DEI.bigWig"

lapply(bw,function(b){
  print(b)
  lapply(c(1:8),function(i){
    print(i)
    res = ALL_lst[[i]]
    lapply(c(1:2),function(n){
      phast <- get_hist(g_rds = split(res, strand(res))[[n]],
                        bw_path = b,
                        seq = 2500,
                        label = paste0(names(ALL_lst)[i],"_",names(split(res, strand(res)))[n])) %>%
        mutate(Histone = gsub("_no_chr.bigWig","",str_split_fixed(b,"/",13)[,13]))
        
      phast
})%>% do.call(rbind,.)->a
  a
  })%>% do.call(rbind,.)->b
  b
})%>% do.call(rbind,.)->tc_lst

tc_lst$type1 <- mapply(function(x) x[1], strsplit(tc_lst$ss, "_"))
tc_lst$type2 <- mapply(function(x) x[2], strsplit(tc_lst$ss, "_"))
tc_lst[tc_lst$type2=="-",]$Pos <- -tc_lst[tc_lst$type2=="-",]$Pos
plot <- tc_lst[,.(mean = mean(y)), by=.(Histone,type1,Pos)]


```
```{r fig.width = 2.4,fig.height=2.5}
c("#FC8D62","#66C2A5","#8DA0CB","#E78AC3" ,"#ED646B","#a9b2b6" ,"#ecd276","#da954f") -> col
names(col) <- c("scATS","SCAFE","CamoTSS","GENCODE","CAGE","Random","TSSr","scTSS")

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/ATAC//all.8.pdf",height = 2.5,width = 2.4)

ggplot(plot, aes(x = Pos, y = mean,colour = factor(type1, levels = names(col)))) + 
    geom_line(size = 0.7,alpha=.7) + 
  scale_color_manual(values =col)+
  ylab("Average signal")+xlab("Distance to TSS")+
    scale_x_continuous(breaks = seq(-2500,2500,1250))+
    guides(shape = FALSE, colour = guide_legend(override.aes = list(size = 3))) + 
    theme_classic()+
    theme(axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          axis.title =  element_text(size = 11),
          legend.text=element_text(size=10),
          plot.title = element_text (hjust = 0.5),
          legend.title = element_blank(),
          strip.text = element_text(size = 14),
          legend.position = c(0.8,0.66))+
  guides(fill = guide_legend( nrow = 5))
dev.off()
```

# g
```{r}
ALL_lst <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/scATSRD_SCAFE_camotss_obj8_0106.Rds")


lapply(c(1:8),function(i){
  print(i)
  res = ALL_lst[[i]]
  lapply(c(1:2),function(n){
      phast <- get_phast(g_rds = split(res, strand(res))[[n]],
                        bw_path = "/mnt/raid66/Personal_data/xuzijie/ref/MusMus/UCSC/mm10.60way.phastCons.bw",
                        seq = 2500,
                        label =paste0(
                          names(ALL_lst)[i],"_",
                          names(split(res,strand(res)))[n]))
      phast
    
  })%>% do.call(rbind,.)->utr_lst
})%>% do.call(rbind,.)->pha_lst

pha_lst$type1 <- mapply(function(x) x[1], strsplit(pha_lst$ss,"_",2))
pha_lst$type2 <- mapply(function(x) x[2], strsplit(pha_lst$ss,"_",2))

pha_lst[pha_lst$type2=="-",]$Pos <- -pha_lst[pha_lst$type2=="-",]$Pos
plot <- pha_lst[,.(mean = mean(y)), by=.(type1,Pos)]

```


```{r fig.width = 2.4,fig.height=2.5}
level <- c("scATS","SCAFE","CamoTSS","TSSr","scTSS","GENCODE","CAGE","Random")

pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/mm_BM/V8/02accuracy/phast/all.8_1125.pdf",height = 2.5,width = 2.4)

ggplot(plot, aes(x = Pos, y = mean,colour = factor(type1,levels = level))) + 
    geom_line(size = 0.7, alpha=.7) +
  scale_color_manual(values = col)+
  ylab("PhastCons elements")+xlab("Distance to TSS")+
    guides(shape = FALSE, colour = guide_legend(override.aes = list(size = 6))) + 
    theme_classic()+
    scale_x_continuous(breaks = seq(-2500,2500,1250))+
    theme(axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          axis.title =  element_text(size = 14),
          legend.text=element_text(size=10),
          plot.title = element_text (hjust = 0.5),
          legend.title = element_blank(),
          strip.text = element_text(size = 14),
          legend.position = c(0.83,0.66)) 

dev.off()
```

