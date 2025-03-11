# Usage
**scATS** can analyze both single-end and paired-end 5'-end scRNA-seq data, enabling direct quantification using Seurat objects, while incorporating **RNA degradation** (RD) modeling through expectation-maximization (EM).

The main functions of scATS are: `TSSCDF`, `FindMarkers`, `psi`, `Sashimi`

## TSS inference and quantification
The input files include:
- Seurat object (R object)
- Alignment file (bam file)
- Annotation file (gtf file)

```r
# load packages
suppressMessages({
library(scATS)
library(Seurat)
})

# load input files
seu <- readRDS("/mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/seurat/gene_wnn_res3_seurat.rds")
Genes <- rownames(seu)
gtfFile <- file.path("/mnt/raid61/Personal_data/tangchao/Document/gencode/mouse/release_M25/gencode.vM25.primary_assembly.annotation.sorted.gtf")
bams <- list.files("/mnt/raid66/Personal_data/xuzijie/task/07ATS/00data/thesis/mm_BM/bam/01bam","bam$",full.names = T)
file.exists(bams)

# quantification using wrapper function
scats <- TSSCDF(object = seu, bam = bams, gtf = gtfFile, genes = Genes, verbose = TRUE)
scats
```

    class: scATSDataSet 
    dim: 11464 388 
    metadata(1): version
    assays(2): counts psi
    rownames(11464): CAAA01118383.1@GL456216.1:16275:+
    Vamp7@GL456233.1:39139:- ...
    Gm47283@chrY:90793346:+ Gm47283@chrY:90793469:+
    rowData names(19): gene_id gene_name ...
    UnclassifiedReads AllReads
    colnames(388): AACCGCGAGAGTCTGG-1
    AACTCAGAGCAGCGTA-1 ... TTTCCTCAGTCAAGGC-1
    TTTCCTCCACCAGATT-1
    colData names(32): orig.ident nCount_RNA ...
    celltype wsnn_cell_type

The quantification results are stored in  `scats@rowRanges` and contain the following columns:  

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| seqnames         | The chromosomal name of the gene to which TSS belongs  |
| ranges           | The genomic locus of TSS, where the growth rate is highest  |
| strand           | The genomic strand of the gene to which TSS belongs   |
| gene_id          | The Ensembl ID of the gene to which TSS belongs   |
| gene_name        | The HGNC Symbol of the gene to which TSS belongs   |
| TSS              | The ID of the TSS in the format of [seqnames:ranges:strand]   |
| Region           | The growth area of the sorted 5'-starts of read 1, or it can also be interpreted as the<br> TSS cluster. Another region must show a significant increase in the number of reads.   |
| Base             | The span length of the start, which is actually included in the calculation, equal to sigma   |
| PSI              | The percent spliced in (PSI) of TSS  |
| Percent          | The ratio of TSS reads   |
| Annotated        | The nearest annotated TSS locus that exists in the region, if it is NA, <br>indicates that there are no annotated TSS loci in the region   |
| Greedy           | Greedy for the first TSS. If it is TRUE, it indicates that the first TSS is quantified<br> even if it does not meet the specified condition."   |
| AUC              | The area under the cumulative distribution curve : close to 1, indicates no degradation,<br>while close to 0 indicates severe degradation.   |
| MaxDist          | The y-axis value corresponding to the farthest point of the cumulative distribution curve:<br> close to 1 indicates no degradation, while close to 0 indicates severe degradation.  |
| AboveRandom      | 在对角线上的read比例   |
| ApicesX          | The x-axis value corresponding to the farthest point of the cumulative distribution curve:<br> close to 0 indicates no degradation, while close to 1 indicates severe degradation.   |
| ApicesY          | Equal to MaxDist.  |
| alpha1           | The degradation index fitted using the EM algorithm: the larger the value, the more severe the degradation.    |
| alpha2           | The degradation index fitted using the EM algorithm, which has now been deprecated.   |
| theta            | The PSI value fitted using the EM algorithm.   |
|UnclassifiedReads | The ratio of discarded reads.   |
| AllReads         | The total number of reads used for quantification.   |

***Note : All the following analyses are based on the `scats` object.***

## Finding differentially expressed ATSs
### Finding markers of all group
```r
DE <- scATS::FindMarkers(object = scats, groupBy = "CellType", group1 = NULL, group2 = NULL, cores = 20, majorOnly = F) 
DE[1:2,]
```
                         TSS  G1    G2 n1  n2 N1  N2 Cells1 Cells2      PSI1
    1 A1BG-AS1@19:58347753:+ AT2 Other 11 270 12 290    357   6826 0.9166667
    2 A1BG-AS1@19:58347753:+ Fib Other  4 277  6 296    385   6798 0.6666667
           PSI2 SudobulkPSI1 SudobulkPSI2  wald.test wilcox.test    prop.test
    1 0.9310345    0.8253968    0.9417671 0.84819606   0.8511922 5.386585e-04
    2 0.9358108    0.6428571    0.9424460 0.02694748   0.0105782 1.232170e-09

### Finding markers of given group
```r
DE <- scATS::FindMarkers(object = scats, groupBy = "CellType", group1 = "AT2", group2 = "Fib", cores = 20, majorOnly = F) 
DE[1:2,]
```
                          TSS     G1     G2    n1    n2    N1    N2 Cells1 Cells2      PSI1
    1: A1BG-AS1@19:58347753:+    AT2    Fib    11     4    12     6    357    385 0.9166667
            PSI2 SudobulkPSI1 SudobulkPSI2 wald.test wilcox.test prop.test
    1: 0.6666667    0.8253968    0.6428571  0.208955   0.2181715 0.1014264

In addition, you can specify the host genes used in the calculation by setting the `gene` parameter.
The `DE` object contains following columns:

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| TSS              | The ID of the TSS in the format of [gene_name:seqnames:ranges:strand].   |
| G1/2             | The levels of group.   |
| n1/2             | The expression of the TSS in the G1/2.   |
| N1/2             | The expression of the host gene in the G1/2.   |
| Cells1/2         | The number of cell-type corresponding to G1/2.   |
| PSI1/2           | The average sum of all individual cell PSIs.   |
| SudobulkPSI1/2   | The PSI value calculated by combining all reads and treating them as a pseudobulk sample.   |
| wald.test        |  psi在两组的差异  |
| wilcox.test      |    |
| prop.test        |    |



## Calculating PSI
```r
psi <- scATS::psi(object = scats, groupBy = "CellType")
psi[1:2,]
```
                         TSS    groupBy Cells   N   n      mean        sd
    1 A1BG-AS1@19:58347753:+          B  2697  67  49 0.7038972 0.4525730
    2 A1BG-AS1@19:58347753:+ Epithilial  4299 142 105 0.7300469 0.4424003
              se         ci median Q1 Q3 mad iqr PseudobulkPSI
    1 0.05529059 0.11039123      1  0  1   0   1     0.7441176
    2 0.03712541 0.07339439      1  0  1   0   1     0.7355705

The `psi` object contains following columns:

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| TSS              | The ID of the TSS in the format of [gene_name:seqnames:ranges:strand].   |
| groupBy          | The level of group.   |
| Cells            | The number of cell-type corresponding to a given groups.   |
| N                | The expression of the host gene in a given groups.   |
| n                | The expression of the TSS in a given groups.   |
| mean             |  都是psi  |
| sd               |    |
| se               |    |
| ci               | 置信区间   |
| median           |    |
| Q1               |    |
| Q3               |    |
| mad              |  std  |
| iqr              |    |
| PseudobulkPSI    | The PSI value calculated by combining all reads and treating them as a pseudobulk sample.   |


## Sashimi plots

```r
# Take SCNN1D gene as an example
scATS::Sashimi(object = scats, 
               bam = bams,
               xlimit = c(1280352, 1282325),
               transcripts = c("SCNN1D-201","SCNN1D-206","SCNN1D-205"),
               gtf = gtfFile, 
               gene = "SCNN1D", 
               TSS = c(1280452, 1281224, 1282125), # TSS位点
               free_y = T,#是否scale
               base_size = 12, #read部分字体大小
               rel_height=0.9, #注释/read ，小于1 read部分比例更大
               fill.color = "grey",
               line.color = "red",
               line.type = 3) -> p
p
```
![sashimi](image.png)

```r
sessionInfo()
```

    R version 3.6.0 (2019-04-26)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 16.04.7 LTS

    Matrix products: default
    BLAS:   /usr/lib/openblas-base/libblas.so.3
    LAPACK: /usr/lib/libopenblasp-r0.2.18.so

    locale:
    [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
    [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    loaded via a namespace (and not attached):
    [1] bitops_1.0-7                matrixStats_1.3.0           bit64_4.0.5                
    [4] RcppAnnoy_0.0.19            RColorBrewer_1.1-3          progress_1.2.3             
    [7] httr_1.4.7                  GenomeInfoDb_1.22.1         tools_3.6.0                
    [10] irlba_2.3.3                 utf8_1.2.4                  R6_2.5.1                   
    [13] KernSmooth_2.23-24          DBI_1.2.3                   BiocGenerics_0.32.0        
    [16] colorspace_2.1-0            tidyselect_1.2.1            prettyunits_1.2.0          
    [19] bit_4.0.5                   curl_5.2.1                  compiler_3.6.0             
    [22] cli_3.6.3                   Biobase_2.46.0              DelayedArray_0.12.3        
    [25] rtracklayer_1.46.0          scales_1.3.0                lmtest_0.9-40              
    [28] ggridges_0.5.6              askpass_1.2.0               rappdirs_0.3.3             
    [31] stringr_1.5.1               digest_0.6.25               Rsamtools_2.2.3            
    [34] R.utils_2.12.3              XVector_0.26.0              pkgconfig_2.0.3            
    [37] htmltools_0.5.0             parallelly_1.37.1           dbplyr_2.5.0               
    [40] fastmap_1.2.0               htmlwidgets_1.5.1           rlang_1.1.1                
    [43] rstudioapi_0.11             RSQLite_2.3.7               generics_0.1.3             
    [46] zoo_1.8-12                  jsonlite_1.8.8              ica_1.0-3                  
    [49] BiocParallel_1.20.1         dplyr_1.1.4                 R.oo_1.26.0                
    [52] RCurl_1.98-1.16             magrittr_2.0.3              GenomeInfoDbData_1.2.2     
    [55] patchwork_1.2.0             Matrix_1.5-3                Rcpp_1.0.13                
    [58] munsell_0.5.1               S4Vectors_0.24.4            fansi_1.0.6                
    [61] reticulate_1.38.0           lifecycle_1.0.4             R.methodsS3_1.8.2          
    [64] stringi_1.8.4               MASS_7.3-57                 SummarizedExperiment_1.16.1
    [67] zlibbioc_1.32.0             Rtsne_0.17                  BiocFileCache_1.10.2       
    [70] grid_3.6.0                  blob_1.2.4                  ggrepel_0.9.5              
    [73] listenv_0.9.1               parallel_3.6.0              crayon_1.5.3               
    [76] lattice_0.20-45             Biostrings_2.54.0           cowplot_1.1.3              
    [79] splines_3.6.0               GenomicFeatures_1.38.2      hms_1.1.3                  
    [82] knitr_1.29                  pillar_1.9.0                igraph_2.0.3               
    [85] GenomicRanges_1.38.0        future.apply_1.11.2         codetools_0.2-20           
    [88] biomaRt_2.42.1              stats4_3.6.0                leiden_0.4.3.1             
    [91] XML_3.99-0.3                glue_1.7.0                  png_0.1-8                  
    [94] vctrs_0.6.5                 gtable_0.3.5                openssl_2.2.0              
    [97] RANN_2.6.1                  future_1.33.2               cachem_1.1.0               
    [100] ggplot2_3.5.1               xfun_0.16                   survival_3.7-0             
    [103] tibble_3.2.1                GenomicAlignments_1.22.1    AnnotationDbi_1.48.0       
    [106] memoise_2.0.1               IRanges_2.20.2              cluster_2.1.6              
    [109] globals_0.16.3              fitdistrplus_1.2-1          ROCR_1.0-11  