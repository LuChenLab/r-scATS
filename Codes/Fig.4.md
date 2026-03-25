---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
suppressMessages({library(Seurat)
library(scATS)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(parallel)
library(RColorBrewer)
library(showtext)})

cellorder2 <- c("B","Mast","Epithilial")
color2 <- c("#aaddcd","#b6c2dd","#f2bddd")
names(color2) <- cellorder2

c("#7AC2E2", "#E56C83") -> co_cancer
order_cancer <- c("Non-malignant","Malignant")
names(co_cancer) <- order_cancer

```

# a
```{r fig.height=3.5, fig.width=7}
cell_de <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.DE_epi.Rds")
cell_de$pad <- p.adjust(cell_de$p)

cell_de$condition <- as.factor(
    ifelse(cell_de$p < 0.05 & abs(cell_de$del_psi) > cutoff,ifelse(cell_de$del_psi > cutoff,'Malignant','Non-malignant'),'Common'))
table(cell_de$condition)

cell_de %>% dplyr::group_by(condition) %>% dplyr::add_count() -> cell_de
cell_de$type2 <- paste0(cell_de$condition, " (", cell_de$n, ")")

cell_de$type2 <- factor(cell_de$type2, levels = c("Non-malignant (1320)", "Common (22371)", "Malignant (2011)"))

saveRDS(cell_de, "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.DE.Rds")

```

```{r fig.height=3.5, fig.width=4.5}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/03infercnv/infercnv.volcano3.pdf", width = 4.5, height = 3.5)

ggplot(cell_de, aes(x = del_psi, y = -log10(p), 
                color = factor(type2, levels = c("Non-malignant (1320)", "Common (22371)", "Malignant (2011)")))) +
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
  geom_vline(xintercept = c(-0.1, 0.1), color = "darkgrey", linetype = "dashed", linewidth = 0.5) +
  ylim(0, 30.5)+ xlim(-0.6,1)+
  geom_point(data = label_data_right,
             aes(x = del_psi, y = -log10(p)),
             color= '#E56C83', size = 4.5, alpha = 0.2) +
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
  xlab("Delt PSI Malignant vs Non-malignant")+
  ylab("-Log10 (P-value)")
  
dev.off()
```

# b
```{r}
library(org.Hs.eg.db)
library(clusterProfiler)
theta_infercnv <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.theta_epi.Rds")
theta_infercnv$gene <- mapply(function(x) x[1], strsplit(theta_infercnv$TSS,"@"))

theta_infercnv[theta_infercnv$type=="Malignant",] -> Malignant_df
Malignant_df$ID <- mapIds(org.Hs.eg.db,
                  keys = as.character(Malignant_df$gene),
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals="first")

ora_res <- enrichGO(gene = Malignant_df$ID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",#这里指定ID类型
                 ont = "BP", # "BP", "MF", "CC" 
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.5,
                 minGSSize = 5,# 最少的基因数量
                 maxGSSize = 300, # 最大的基因数量
                 readable = T # 把ENTREZID转换为SYMBOL
                 ) # 用什么函数富集都可以，主要是整理成 ora_res2 的格式

saveRDS(ora_res, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.theta_go.Rds")
```

```{r}
library(aPEAR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(cols4all)


ora_res <- readRDS(file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.theta_go.Rds")
ora_res2=clusterProfiler::simplify(ora_res,
                                         cutoff=0.5,
                                         by="p.adjust",
                                         select_fun=min)
plot <- ora_res@result[order(ora_res@result$pvalue, decreasing = F),]
plot2 <- plot[plot$ID %in% c(
                            plot[grep("kappaB",plot$Description),]$ID,
                            plot[grep("T cell",plot$Description),]$ID[1:20],
                             plot[grep("epithelial",plot$Description),]$ID[1:20],
                             plot[grep("lung",plot$Description),]$ID,
                             plot[grep("virus",plot$Description),]$ID),]
plot2$pvalue2 <- -log10(plot2$pvalue)
```


```{r fig.height=5, fig.width=8}
library(RColorBrewer)
library(ggplot2)
set.seed(1)
p <- enrichmentNetwork(plot2,
                        colorBy = 'pvalue',
                        colorType = c("pval"),
                        nodeSize = 'Count',
                        fontSize = 4,
                        drawEllipses = T,
                        outerCutoff=0.3,
                       minClusterSize=3,
                        pCutoff = -50,
                        verbose = TRUE) +
  scale_color_gradient(low = "lightgrey", high = "red")
p


```


```{r fig.height=5, fig.width=8}
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/02DE/NSCLC.infercnv.theta_go.pdf", height = 5, width = 8)
p
dev.off()
```


# d
## LRS models
```{r}

data_df <- readRDS( "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/feature_v5.Rds") 
sizeN = 129

data_df <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/feature_v5.Rds")
data_df[data_df$Label==1,]

## 确定负集
data_df <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/feature_v5.Rds")
sizeN = 129
set.seed(2250)
# rownames(data_df) <- data_df$TSS
# data_df$TSS <- NULL

data <- data_df[,-ncol(data_df)]
groudTruth <- data[rownames(data) %in% rownames(data_df[data_df$Label==1,]),]
Unlabeled <-  data[rownames(data) %in% rownames(data_df[data_df$Label==0,]),]

##PU-learning Spy
indx <- sample(2, nrow(groudTruth), replace = T, prob = c(0.9, 0.1))
s <- groudTruth[indx == 2,] # 10% for spy
Ps <- groudTruth[indx == 1,]
Us <- rbind(Unlabeled, s )

Ps$NewLabel <- 1
Us$NewLabel <- 0

### oversample
rn_train <- rbind(Ps, Us)
balance.over <- ovun.sample(NewLabel~., data = rn_train, p = 0.5, seed = 1, method = "over")$data #-40

table(balance.over$NewLabel)

## naiveBayes classify

classifier <- naiveBayes(NewLabel ~., balance.over)
rn_pred_u <- predict(classifier, Unlabeled, type = "raw" )
rn_pred_s <- predict(classifier, s, type = "raw" )

set.seed(2250)
tr <- quantile(rn_pred_s[,2], probs = seq(0,1,0.01))[2][[1]] # 1%概率为正样本的值
index_RN <- which(rn_pred_u[,2] < tr) # 小于这个值便为负样本
RN <- Unlabeled[sample(index_RN,sizeN,replace = F),]
RN$Label <- 0
PN <- groudTruth %>% dplyr::mutate(Label = 1)

tr;length(index_RN)
```

```{r}
### select reliable negative

data_df$Label <- ifelse(rownames(data_df) %in% rownames(PN), "Pos",
                        ifelse(rownames(data_df) %in% rownames(RN), "Neg", "Unlabel"))

save(RN, PN, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/PN_RN.RData")
```


```{r}
## model training
load("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/PN_RN.RData")

set.seed(2250)
Dataset <- rbind(RN,PN)

Dataset$Label <- as.factor( Dataset$Label)
folds <- createFolds(y=Dataset$Label,k=10)  ### 10-fold
RF_classifier <- list()
roc_obj <- list()
auc_value <- list()


for(i in 1:10){
  fold_test <- Dataset[folds[[i]],] #folds[[i]] as test data
  fold_train <- Dataset[-folds[[i]],] # training data
  
  RF_classifier_v1 <- randomForest(Label ~ ., data = fold_train, mtry = 5, ntree = 1000, importance = T)
  RF_classifier[[i]] <- RF_classifier_v1
  
  ##### validation and ROC
  pred_rf_v1 <- predict(RF_classifier_v1, newdata = fold_test[,-ncol(fold_test)], type = "prob")
  obs_p_rf_v1 <- data.frame(pred = pred_rf_v1[,2],
                            obs = fold_test$Label )
  rf_roc_v1 <- roc(obs ~ pred, obs_p_rf_v1, levels = c("0", "1"))
  auc_value_v1 <- as.numeric(auc(rf_roc_v1))
  
  roc_obj[[i]] <- rf_roc_v1
  auc_value[[i]] <- auc_value_v1
}
```

```{r fig.height=3,fig.width=3.5}
auc_mean <- c()
for (i in 1:10){
  auc_mean <- append(auc_mean, auc(roc_obj[[i]]))
}
mean(auc_mean)

AUC <- ggroc(roc_obj, alpha = 0.5, linetype = 1, size = 0.8, legacy.axes = TRUE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 14)) +
  annotate('text', x=0.7, y=0.3, size=5, color='red', label = paste0("AUCmean = ",round(mean(auc_mean),2))) +
  labs(x = "1-Specificity", y = "Sensitivity")

cairo_pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5//Fig.AUC_10fold.pdf", width = 3.5,height = 3)
AUC
dev.off()

save(RF_classifier, roc_obj, auc_value, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/model_v1.Rdata")
```


```{r}
load("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/PN_RN.RData")

Dataset <- rbind(RN,PN)
write.csv(Dataset, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/Pos_Neg.csv", col.names = F, row.names = F, quote = F)
write.csv(data_df[,-ncol(data_df)], file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/all_ATS.csv", col.names = F, row.names = F, quote = F)

write.csv(data_df[,-ncol(data_df)], file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/all_ATS_github.csv", col.names = F, row.names = T, quote = F)

```


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


# g
```{r}
df <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/RF.score.Rds")
library(clusterProfiler)
library(org.Hs.eg.db)


df[df$GeneID %in% ats,c("GeneID","pred_label")]
mydf <- df[,c("GeneID","pred_label")]
mydf$Gene <- stringr::str_split_fixed(as.character(mydf$GeneID), "@",2)[,1]
mydf <- mydf[mydf$pred_label=="Predicted functional",]

mydf$ID <- mapIds(org.Hs.eg.db,
                  keys = mydf$Gene,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals="first")

go_pos <- enrichGO(gene = mydf$ID,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",#这里指定ID类型
                   ont = "BP", # "BP", "MF", "CC" 
                   pvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.1,
                   minGSSize = 5,
                   readable = T # 把ENTREZID转换为SYMBOL
                   )

saveRDS(go_pos, file = "/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/RF_ablation_Pos.GO.Rds")

```


```{r}
go_pos <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/RF_ablation_Pos.GO.Rds")
library(gground)
library(ggprism)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

ego_readable <- setReadable(go_pos, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
GO <- as.data.frame(ego_readable)

use_pathway <- group_by(GO, ONTOLOGY) %>%
  top_n(10, wt = pvalue) %>%
  group_by(pvalue) %>%
  top_n(1, wt = Count) %>%
  
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF')))) %>%
  dplyr::arrange(ONTOLOGY, pvalue) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')


width <- 0.5
# x 轴长度
xaxis_max <- max(-log10(use_pathway$pvalue)) + 1
# 左侧分类标签数据
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )
# 绘制富集通路图
plot_enrichment <- function() {
  p <- use_pathway %>%
    ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), width =0.6, alpha =0.8
    ) +
    geom_text(
      aes(x =0.05, label = Description),
      hjust =0, size =5
    ) +
    geom_text(
      aes(x =0.1, label = substr(geneID,start = 1, stop = 50), colour = ONTOLOGY),
      hjust =0, vjust =2.6, size =3.5, fontface ='italic',
      show.legend =FALSE
    ) +
    # 基因数量
    geom_point(
      aes(x = -width, size = Count),
      shape =21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name ='Count', range = c(5,12)) +
    # 分类标签
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2,'mm'),
      inherit.aes =FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes =FALSE
    ) +
    geom_segment(
      aes(x =0, y =0, xend = xaxis_max, yend =0),
      linewidth =1.5,
      inherit.aes =FALSE
    ) +
    labs(y =NULL) +
    scale_fill_manual(name ='Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max,2),
      expand = expansion(c(0,0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
return(p)
}
```


```{r fig.height=7, fig.width=6}

pal <- c('#7bc4e2', '#acd372', '#fbb05b', '#ed6ca4')
pal <- c('#eaa052', '#b74147', '#90ad5b', '#23929c')
pal <- c('#c3e1e6', '#f3dfb7', '#dcc6dc', '#96c38e')

enrichment_plot <- plot_enrichment()
pdf("/mnt/raid66/Personal_data/xuzijie/task/07ATS/02result/thesis/NSCLC/V9/05model/xzj/v5/RF_ablation_Pos.go.pdf", width = 6, height = 7)

enrichment_plot + ggtitle("GO enrichment of 9,035 positive TSS host genes")
dev.off()

```