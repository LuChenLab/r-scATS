
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

# e

```{r}
features <- read.csv("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/RF.feature.evaluate_ablation1.csv", stringsAsFactors = F)
range(features$AUROC)
range(features$AUC.PR)

```

```{r fig.height=3, fig.width=5}
features$mean <- (features$AUROC+features$AUC.PR) /2
colnames(features) <- c("AUROC","AUC-PR","num","Mean")
lapply(c(1,2,4), function(i){
  ggplot(features,
         aes(x= num, 
          y=features[,i])) +
  geom_line(color = "#cccccc")+
    geom_point(color = "#87CEEB") +
  theme_bw() + labs(x="Number of deleted features",y=colnames(features)[i], title = "")+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) + ylim(0.4,1) + xlim(0,33)+
    geom_vline(xintercept = c(16), color = "red", linetype = "dotted")
}) -> p_lst
which(diff(features$AUROC)>0)
```


```{r fig.height=3, fig.width=14}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/RF.feature.dele.pdf", width = 14, height = 3)

cowplot::plot_grid(p_lst[[1]], p_lst[[2]],p_lst[[3]],  nrow=1)

dev.off()
```
