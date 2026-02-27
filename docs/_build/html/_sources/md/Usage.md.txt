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
seu <- readRDS("demo_seurat.Rds")
Genes <- rownames(seu)
gtfFile <- file.path("demo.gtf")
bams <- "demo.bam"
file.exists(bams)

# quantification using wrapper function
scats <- TSSCDF(object = seu, bam = bams, gtf = gtfFile, genes = Genes, verbose = TRUE)
scats

```

    class: scATSDataSet
    dim: 2043 2000
    metadata(2): version parameters
    assays(3): counts psi theta
    rownames(2043): OR4F5@1:69063:+ OR4F5@1:69071:+ ... OR2G6@1:248508175:+
    ZNF672@1:248838224:+
    rowData names(12): TSS gene_id ... alpha theta
    colnames(2000): TGGACGCTCCTTCAAT-1 CGAGCACGTCAGAATA-1 ...
    ATCGAGTAGCGGCTTC-1 CTAGAGTAGTGCCATT-1
    colData names(4): orig.ident nCount_RNA nFeature_RNA Cell



The quantification results are stored in  `scats@rowRanges` and contain the following columns:  

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| seqnames         | The chromosomal name of the gene to which TSS belongs.  |
| ranges           | The genomic locus of TSS, where the growth rate is highest.  |
| strand           | The genomic strand of the gene to which TSS belongs.   |
| gene_id          | The Ensembl ID of the gene to which TSS belongs.   |
| gene_name        | The HGNC Symbol of the gene to which TSS belongs.   |
| TSS              | The ID of the TSS in the format of [seqnames:ranges:strand].   |
| Region           | The growth area of the sorted 5'-starts of read 1, or it can also be interpreted as the<br> TSS cluster. Another region must show a significant increase in the number of reads.   |
| PSI              | The percent spliced in (PSI) of TSS.  |
| Percent          | The ratio of TSS reads.   |
| Annotated        | The nearest annotated TSS locus that exists in the region, if it is NA, <br>indicates that there are no annotated TSS loci in the region.   |
| Greedy           | Greedy for the first TSS. If it is TRUE, it indicates that the first TSS is quantified<br> even if it does not meet the specified condition."   |
| beta             | The area under the cumulative distribution curve : close to 1, indicates no degradation,<br>while close to 0 indicates severe degradation.   |
| alpha            | The degradation index fitted using the EM algorithm: the larger the value, the more severe the degradation.    |
| theta            | The PSI value fitted using the EM algorithm.   |
| AllReads         | The total number of reads used for quantification.   |

***Note : All the following analyses are based on the `scats` object.***


## Finding differentially expressed ATSs
### Finding markers of all group

**Identifying ATS markers based on `θ` value.**

```r
DE <- scATS::FindMarkersByTheta(object = scats, groupBy = "CellType", group1 = NULL, group2 = NULL, cores = 20)
DE[1:2,]
```

          gene            TSS group1    theta1    theta2
        <char>         <char> <char>     <num>     <num>
    1:  OR4F5 OR4F5@1:69063:+      A 0.2780345 0.2835192
    2:  OR4F5 OR4F5@1:69063:+      B 0.2862817 0.2620500
        cell1    cell2   percent1   percent2          p
        <int>    <int>      <num>      <num>      <num>
    1:      6       17   0.896861  1.2772352  0.9682337
    2:     10       13   1.459854  0.9885932  0.5724791

The `DE` object contains following columns:

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| gene             | The HGNC Symbol of the gene to which TSS belongs.   |
| TSS              | The ID of the TSS in the format of [gene_name:seqnames:ranges:strand].   |
| group1           | The levels of group.   |
| theta1/2         | The theta value of the TSS in the group.   |
| cell1/2          | The number of cell-type corresponding to group.   |
| percent1/2       | The ratio of TSS reads.  |
| p                | The p-value from the Wilcoxon test for differences in psi values among groups.   |


**Identifying ATS markers based on `ψ` value.**

```r
DE <- scATS::FindMarkers(object = scats, groupBy = "CellType", group1 = NULL, group2 = NULL, cores = 20, majorOnly = F) 
DE[1:2,]
```
                   TSS     G1     G2    n1    n2    N1    N2 Cells1 Cells2      PSI1
                <char> <char> <char> <int> <int> <int> <int>  <int>  <int>     <num>
    1: OR4F5@1:69071:+      B  Other    41    83   156   277    685   1315 0.2088708
    2: OR4F5@1:69156:+      B  Other    40    90   156   277    685   1315 0.2169872
            PSI2 PseudobulkPSI1 PseudobulkPSI2  wald.test wilcox.test    prop.test
            <num>          <num>          <num>      <num>       <num>        <num>
    1: 0.2488509      0.3492958      0.2088773 0.42343003  0.32714124 3.331726e-11
    2: 0.2964099      0.1309859      0.2663185 0.06632703  0.08046456 7.544541e-12

The `DE` object contains following columns:

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| TSS              | The ID of the TSS in the format of [gene_name:seqnames:ranges:strand].   |
| G1/2             | The levels of group.   |
| n1/2             | The expression number of the TSS in the G1/2.   |
| N1/2             | The expression number of the host gene in the G1/2.   |
| Cells1/2         | The number of cell-type corresponding to G1/2.   |
| PSI1/2           | The average sum of all individual cell PSIs.   |
| PseudobulkPSI1/2   | The PSI value calculated by combining all reads and treating them as a pseudobulk sample.   |
| wald.test        | The p-value from the Wald test for differences in psi values among groups.  |
| wilcox.test      | The p-value from the Wilcoxon test for differences in psi values among groups.   |
| prop.test        | The p-value from the Proportion test for differences in psi values among groups.   |



### Finding markers of given group

```r
### based on θ value
DE <- scATS::FindMarkersByTheta(object = scats, groupBy = "CellType", group1 = NULL, group2 = NULL, cores = 20)
### based on  ψ value
DE <- scATS::FindMarkers(object = scats, groupBy = "CellType", group1 = "A", group2 = "B", cores = 20, majorOnly = F,gene = "OR4F5") 
```

In addition, you can specify the host genes used in the calculation by setting the `gene` parameter.


## Calculating PSI

**Calculating PSI based on `θ` value.**
```r
theta <- scATS::ThetaByGroup(object = scats, gene = "OR4F5",groupBy = "CellType")
theta[1:2,]
```

    Group       TSS     alpha     theta
    <char>    <char>     <num>     <num>
    1:      A 1:69156:+ 0.1061093 0.1718889
    2:      A 1:69090:+ 0.3006347 0.2107304


**Calculating PSI based on `ψ` value.**

```r
psi <- scATS::psi(object = scats, groupBy = "CellType", TSS=c("OR4F5@1:69071:+", "OR4F5@1:69156:+"))
psi[1:2,]
```
                    TSS groupBy Cells   N  n    mean        sd
    1 OR4F5@1:69071:+       A   669 129 39 0.2562320 0.4188128
    2 OR4F5@1:69071:+       B   685 156 41 0.2088708 0.3822363
              se         ci median Q1  Q3 mad iqr PseudobulkPSI
    1 0.03687441 0.07296233      0  0 0.5   0 0.5     0.2067138
    2 0.03060340 0.06045356      0  0 0.2   0 0.2     0.3492958


The `theta` or `psi` object contains following columns:

| column name      |                                           content                                          |
| ---------------- | ------------------------------------------------------------------------------------------ |
| Group/groupBy    | The level of group.   |
| TSS              | The ID of the TSS in the format of [gene_name:seqnames:ranges:strand].   |
| alpha            | The degradation index fitted using the EM algorithm.   |
| Cells            | The number of cell-type corresponding to a given groups (The same applies to the following.).   |
| N                | The expression of the host gene.   |
| n                | The expression of the TSS.   |
| mean             | The mean of PSI.  |
| sd               | The standard deviation (std) of PSI.  |
| se               | The standard error (SE) of PSI.   |
| ci               | The confidence interval (CI) of PSI.   |
| median           | The median of PSI.   |
| Q1               | The first quartile (Q1) of PSI.   |
| Q3               | The third quartile (Q3) of PSI.   |
| mad              | The median absolute deviation (MAD) of PSI.  |
| iqr              | The interquartile range (IQR) of PSI.   |
| theta/PseudobulkPSI    | The PSI value calculated by combining all reads and treating them as a pseudobulk sample.   |


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

```r
    R version 4.5.1 (2025-06-13)
    Platform: x86_64-conda-linux-gnu
    Running under: CentOS Linux 7 (Core)

    Matrix products: default
    BLAS/LAPACK: /mnt/data4/xuzijie/miniconda3/envs/r451/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

    locale:
    [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C               LC_TIME=zh_CN.UTF-8
    [4] LC_COLLATE=zh_CN.UTF-8     LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=zh_CN.UTF-8
    [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                  LC_ADDRESS=C
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C

    time zone: America/New_York
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] scATS_0.5.5

    loaded via a namespace (and not attached):
    [1] RcppAnnoy_0.0.19            splines_4.5.1               later_1.3.2
    [4] bitops_1.0-7                tibble_3.2.1                R.oo_1.24.0
    [7] polyclip_1.10-0             graph_1.68.0                xts_0.13.1
    [10] XML_3.99-0.14               rpart_4.1-15                lifecycle_1.0.4
    [13] OrganismDbi_1.32.0          ensembldb_2.14.1            globals_0.16.3
    [16] lattice_0.22-9              MASS_7.3-58.3               backports_1.4.1
    [19] magrittr_2.0.3              Hmisc_4.6-0                 plotly_4.10.4
    [22] httpuv_1.6.15               Seurat_4.0.5                sctransform_0.3.2
    [25] spatstat.core_2.3-2         askpass_1.2.0               spatstat.sparse_3.0-3
    [28] reticulate_1.24             ggbio_1.38.0                cowplot_1.1.3
    [31] pbapply_1.7-2               DBI_1.2.2                   RColorBrewer_1.1-3
    [34] abind_1.4-5                 zlibbioc_1.36.0             Rtsne_0.17
    [37] GenomicRanges_1.42.0        mixtools_1.2.0              purrr_1.0.2
    [40] R.utils_2.11.0              AnnotationFilter_1.14.0     biovizBase_1.38.0
    [43] BiocGenerics_0.40.0         RCurl_1.98-1.5              nnet_7.3-16
    [46] VariantAnnotation_1.36.0    rappdirs_0.3.3              GenomeInfoDbData_1.2.4
    [49] IRanges_2.24.1              S4Vectors_0.28.1            ggrepel_0.9.5
    [52] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-3
    [55] goftest_1.2-3               fitdistrplus_1.1-6          parallelly_1.37.1
    [58] smoother_1.1                leiden_0.3.9                codetools_0.2-19
    [61] DelayedArray_0.16.3         xml2_1.3.6                  tidyselect_1.2.1
    [64] matrixStats_1.2.0           stats4_4.5.1                BiocFileCache_1.14.0
    [67] base64enc_0.1-3             GenomicAlignments_1.26.0    jsonlite_1.8.8
    [70] Formula_1.2-4               ggridges_0.5.3              survival_3.5-5
    [73] segmented_1.3-4             tools_4.5.1                 progress_1.2.3
    [76] ica_1.0-2                   Rcpp_1.0.12                 glue_1.7.0
    [79] gridExtra_2.3               xfun_0.43                   mgcv_1.8-42
    [82] TTR_0.24.3                  MatrixGenerics_1.2.1        GenomeInfoDb_1.30.1
    [85] dplyr_1.1.4                 BiocManager_1.30.22         fastmap_1.1.1
    [88] GGally_2.1.2                latticeExtra_0.6-29         fansi_1.0.6
    [91] openssl_2.1.1               digest_0.6.35               R6_2.5.1
    [94] mime_0.12                   colorspace_2.1-0            scattermore_1.2
    [97] tensor_1.5                  dichromat_2.0-0             jpeg_0.1-9
    [100] spatstat.data_3.0-4         biomaRt_2.46.3              RSQLite_2.2.9
    [103] R.methodsS3_1.8.1           utf8_1.2.4                  tidyr_1.3.1
    [106] generics_0.1.3              data.table_1.15.4           rtracklayer_1.50.0
    [109] prettyunits_1.2.0           httr_1.4.7                  htmlwidgets_1.6.4
    [112] uwot_0.1.16                 pkgconfig_2.0.3             gtable_0.3.5
    [115] blob_1.2.4                  lmtest_0.9-39               XVector_0.30.0
    [118] htmltools_0.5.8.1           RBGL_1.66.0                 ProtGenerics_1.22.0
    [121] SeuratObject_4.0.4          scales_1.3.0                Biobase_2.50.0
    [124] png_0.1-8                   rstudioapi_0.16.0           knitr_1.46
    [127] reshape2_1.4.4              checkmate_2.0.0             nlme_3.1-162
    [130] curl_5.2.1                  cachem_1.0.8                zoo_1.8-9
    [133] stringr_1.5.1               KernSmooth_2.23-20          parallel_4.5.1
    [136] miniUI_0.1.1.1              foreign_0.8-81              AnnotationDbi_1.52.0
    [139] pillar_1.9.0                grid_4.5.1                  reshape_0.8.8
    [142] vctrs_0.6.5                 RANN_2.6.1                  VGAM_1.1-5
    [145] promises_1.3.0              dbplyr_2.5.0                xtable_1.8-4
    [148] cluster_2.1.4               htmlTable_2.3.0             GenomicFeatures_1.42.3
    [151] cli_3.6.3                   compiler_4.5.1              Rsamtools_2.6.0
    [154] rlang_1.1.4                 crayon_1.5.2                future.apply_1.11.2
    [157] mclust_6.0.0                plyr_1.8.9                  stringi_1.8.3
    [160] viridisLite_0.4.2           deldir_1.0-6                BiocParallel_1.24.1
    [163] munsell_0.5.1               Biostrings_2.58.0           lazyeval_0.2.2
    [166] spatstat.geom_3.2-9         Matrix_1.5-3                BSgenome_1.58.0
    [169] hms_1.1.3                   patchwork_1.2.0             bit64_4.0.5
    [172] future_1.33.2               ggplot2_3.5.1               shiny_1.8.1.1
    [175] SummarizedExperiment_1.20.0 kernlab_0.9-29              ROCR_1.0-11
    [178] igraph_1.5.1                memoise_2.0.1               bit_4.0.5
```