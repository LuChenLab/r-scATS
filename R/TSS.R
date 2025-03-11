#' @name TSSCDF
#' @title TSS inference and quantification
#' @importFrom rtracklayer readGFFAsGRanges
#' @importFrom S4Vectors mcols
#' @importFrom pbapply pblapply
#' @importFrom parallel mclapply
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicFeatures makeTxDbFromGFF fiveUTRsByTranscript
#' @import GenomicRanges
#' @import SummarizedExperiment
#'
#' @param object A Seurat object.
#' @param bam,sep Path of bam file(s), and the character string to separate the bams.
#' @param genes Genes for ATS inference and quantification.
#' @param gtf Path of GTF file or GRanges file of GTF file.
#' @param txdb A \code{TxDb} object.
#' @param UTROnly Only infer ATS within annotated 5'-UTR regions or first exons of transcripts without 5'-UTR.
#' @param MinGeneReads The minimum reads of a gene for TSS inference.
#' @param mapqFilter,isSecondaryAlignment,isSupplementaryAlignment,isDuplicate parameters for \code{ScanBamParam} to read bam file(s).
#' @param min.TSS.percent The minimum percent of a TSS.
#' @param min.local.percent The minimum difference between observed and local average derivative of CDF of a TSS.
#' @param p.cutoff The cutoff of p-value of a TSS.
#' @param window The window for local mean value estimating.
#' @param greedy greedy for the first TSS? If this parameter is set to TRUE, the first TSS will be printed even if it does not satisfy the specified condition.
#' @param cores The number of cores for parallel working.
#' @param verbose A logical controlling if a text progress bar is displayed.
#' @param project.name The name of project.
#' @export

TSSCDF <- function (object, bam, sep = "_", genes, gtf, txdb = NULL, UTROnly = FALSE, MinGeneReads = 500,
                    mapqFilter = 255, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = NA,
                    min.TSS.percent = 5, min.local.percent = 1, p.cutoff = 0.01, window = 10, greedy = TRUE,
                    cores = 1L, verbose = FALSE, project.name = NULL)
{
  stopifnot(is(object, "Seurat") | is(object, "SimpleSeurat"))
  if (is(object, "Seurat")) object <- SimplifySeurat(object)
  stopifnot(all(file.exists(bam)))
  stopifnot(all(file.exists(paste0(bam, ".bai"))))
  stopifnot(identical(slot(object, "Cell"), row.names(slot(object, "meta.data"))))

  stopifnot(min.TSS.percent > 0 & min.TSS.percent < 100)
  stopifnot(min.local.percent > 0 & min.local.percent < 100)
  stopifnot(p.cutoff > 0 & p.cutoff < 1)
  stopifnot(window > 0 & window < 100)
  stopifnot(is.logical(greedy))

  if (!file.exists(gtf)) {
    stop("GTF file is not exists.")
  } else {
    if (is.null(txdb)) {
      txdb <- suppressMessages(suppressWarnings(GenomicFeatures::makeTxDbFromGFF(gtf)))
    } else {
      if (!is(txdb, "TxDb")) {
        txdb <- suppressMessages(suppressWarnings(GenomicFeatures::makeTxDbFromGFF(gtf)))
      }
    }
    seqs <- unique(do.call(c, lapply(bam, function(x) names(Rsamtools::scanBamHeader(x)[[1]]$target))))
    GenomicFeatures::isActiveSeq(txdb) <- GenomeInfoDb::seqlevels(txdb) %in% seqs
    gtf <- rtracklayer::readGFFAsGRanges(gtf)
    gtf <- keepSeqlevels(gtf, intersect(GenomeInfoDb::seqlevels(gtf), seqs), pruning.mode = "coarse")
    fiveUTR <- unlist(GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE))
    firstEx <- subset(unlist(GenomicFeatures::exonsBy(txdb, "tx", use.names = TRUE)), exon_rank == 1)
    firstEx <- firstEx[!names(firstEx) %in% names(fiveUTR)]
    fiveUTR <- c(fiveUTR, firstEx)
    fiveUTR$gene_name <- plyr::mapvalues(names(fiveUTR), gtf$transcript_id, gtf$gene_name, warn_missing = F)
  }
  EBG <- GenomicFeatures::exonsBy(txdb, by = "gene")
  IBG <- unlist(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE))
  names(IBG) <- plyr::mapvalues(names(IBG), with(as.data.frame(gtf), transcript_id[type == "transcript"]), with(as.data.frame(gtf), gene_id[type == "transcript"]), warn_missing = FALSE)
  IBG <- unique(IBG)
  stopifnot(all(c("gene_name", "gene_id", "type") %in% names(S4Vectors::mcols(gtf))))
  if (length(bam) > 1) {
    stopifnot(!is.null(sep))
    stopifnot(!is.na(sep))
    SampleName <- sort(gsub(".bam", "", basename(bam)))
    stopifnot(identical(SampleName, sort(unique(substr(slot(object, "Cell"), 20, max(nchar(slot(object, "Cell"))))))))
  }
  if(!is.null(UTROnly)) stopifnot(is.logical(UTROnly))
  stopifnot(is.logical(verbose))
  stopifnot(MinGeneReads > 0)
  stopifnot(mapqFilter >= 0)
  stopifnot(is.integer(as.integer(cores)))
  if (is.null(genes)) genes <- slot(object, "Gene")
  genes <- genes[genes %in% unique(gtf$gene_name)]
  if (length(genes) == 0) stop("There is no gene for ATS identification!")
  Tars <- gtf[with(data.frame(S4Vectors::mcols(gtf)), type == "gene" & gene_name %in% genes)]
  Tars <- Tars[!with(data.frame(S4Vectors::mcols(Tars)), gene_name %in% names(which(table(Tars$gene_name) > 1)))]
  if (length(Tars) == 0) stop("There is no gene for ATS identification!")
  # save.image(file = "/mnt/raid61/Personal_data/tangchao/Temp/image.RData")
  if (verbose) {
    ATS_count <- pbapply::pblapply(X = seq_along(Tars),
                                   function(i) {
                                     tryCatch(gmmc4(object = object, region = Tars[i],
                                                    bams = bam, gtfrange = gtf, UTR = fiveUTR, EBG = EBG, IBG = IBG,
                                                    UTROnly = UTROnly, MinGeneReads = MinGeneReads,
                                                    min.TSS.percent = min.TSS.percent, min.local.percent = min.local.percent,
                                                    p.cutoff = p.cutoff, window = window, greedy = greedy,
                                                    mapqFilter = mapqFilter, isSecondaryAlignment = isSecondaryAlignment,
                                                    isSupplementaryAlignment = isSupplementaryAlignment,
                                                    isDuplicate = isDuplicate, sep = sep),
                                              error = function(e) e)
                                   }, cl = cores)
  } else {
    ATS_count <- parallel::mclapply(X = seq_along(Tars),
                                    function(i) {
                                      tryCatch(gmmc4(object = object, region = Tars[i],
                                                     bams = bam, gtfrange = gtf, UTR = fiveUTR, EBG = EBG, IBG = IBG,
                                                     UTROnly = UTROnly, MinGeneReads = MinGeneReads,
                                                     min.TSS.percent = min.TSS.percent, min.local.percent = min.local.percent,
                                                     p.cutoff = p.cutoff, window = window, greedy = greedy,
                                                     mapqFilter = mapqFilter, isSecondaryAlignment = isSecondaryAlignment,
                                                     isSupplementaryAlignment = isSupplementaryAlignment,
                                                     isDuplicate = isDuplicate, sep = sep),
                                               error = function(e) e)
                                    }, mc.cores = cores)
  }
  ATS_count <- ATS_count[!mapply(is.null, ATS_count)]
  res_class <- mapply(function(x) class(x)[1], ATS_count)
  ErrorMessage <- ATS_count[grepl("[Ee]rror", res_class)]
  if (length(ErrorMessage) > 0) {
    ErrorMessage <- mapply(ErrorMessage, FUN = function(x) if (is(x, "simpleError")) unlist(strsplit(x$message, ":"))[1] else x)
    print(as.data.frame(table(ErrorMessage)))
  }
  ATS_count <- ATS_count[!grepl("[Ee]rror", res_class)]
  if (length(ATS_count) == 0) {
    stop("There is no gene after filtering!")
  }
  paras <- data.frame(do.call(rbind, lapply(ATS_count, function(x) x$paras)))
  row.names(paras) <- with(paras, paste(gene_name, TSS, sep = "@"))
  count <- as.matrix(do.call(rbind, lapply(ATS_count, function(x) x$mat_count)))
  psi <- as.matrix(do.call(rbind, lapply(ATS_count, function(x) x$mat_psi)))
  rR <- as(with(paras, TSS), "GRanges")
  S4Vectors::mcols(rR) <- paras
  names(rR) <- with(paras, paste(gene_name, TSS, sep = "@"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = count, psi = psi), rowRanges = rR, colData = slot(object, "meta.data"))
  scATSDataSet(se = se, reductions = slot(object, "reductions"), project.name = project.name)
}


#' @title FindMarkers
#' @name FindMarkers
#' @description Find markers of given group
#' @importFrom SummarizedExperiment colData
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment rowData
#' @importFrom data.table as.data.table
#' @importFrom parallel mclapply
#'
#' @param object An scATSDataSet from TSS.
#' @param gene genes for FindMarkers.
#' @param groupBy The column name of group in scATSDataSet object meta.data.
#' @param group1,group2 Only for given group.
#' @param ratio1,ratio2 Cut off of cells expressed given TSS.
#' @param deltaPSI The minimum delta PSI between compared groups.
#' @param majorOnly Whether only the predominant TSS over all cells of each gene is considered, or not.
#' @param method Method for statistical test. Currently, only wald.test, prop.test and wilcox.test are supported.
#'
#' @export

FindMarkers <- function(object, gene = NULL, groupBy, group1 = NULL, group2 = NULL, ratio1 = 0.1, ratio2 = 0.1, deltaPSI = 0.1, majorOnly = TRUE, method = NULL, cores = 1) {
  stopifnot(is(object, "scATSDataSet"))
  stopifnot(groupBy %in% colnames(SummarizedExperiment::colData(object)))
  if(!is.null(group1)) stopifnot(length(group1) == 1)
  if(!is.null(group2)) stopifnot(length(group2) == 1)
  stopifnot(ratio1 >= 0)
  stopifnot(ratio1 <= 1)
  stopifnot(ratio2 >= 0)
  stopifnot(ratio2 <= 1)
  stopifnot(deltaPSI >= 0)
  stopifnot(deltaPSI <= 1)
  stopifnot(is.logical(majorOnly))
  if(!is.null(method)) {
    stopifnot(length(method) == 1)
    stopifnot(is.element(method, c("wilcox.test", "prop.test", "wald.test")))
  }

  if(!is.null(group1) & !is.null(group2)) {
    stopifnot(all(c(group1, group2) %in% SummarizedExperiment::colData(object)[, groupBy]))
    object <- object[, SummarizedExperiment::colData(object)[, groupBy] %in% c(group1, group2)]
  }

  gs <- unique(SummarizedExperiment::colData(object)[, groupBy])
  stopifnot(length(gs) > 1)
  stopifnot(all(c(group1, group2) %in% gs))

  if(is.null(group1)) {
    stopifnot(any(gs != "G0"))
    groups <- lapply(gs, function(i) {
      gs[gs != i] <- "G0"
      return(gs)
    })
  } else {
    if(is.null(group2)) {
      stopifnot(is.element(group1, gs))
      groups <- lapply(setdiff(gs, group1), function(i) {
        gs[!gs %in% c(group1, i)] <- NA
        return(gs)
      })
    } else {
      groups <- gs
      groups[!groups %in% c(group1, group2)] <- NA
      groups <- list(groups)
    }
  }

  if(is.null(gene)) {
    # if(is.null(TSSs)) TSSs <- row.names(object)
    TSSs <- row.names(object)
  } else {
    TSSs <- row.names(object[with(as.data.frame(SummarizedExperiment::rowData(object)), gene_name %in% gene)])
  }
  TSSs <- TSSs[TSSs %in% row.names(object)]
  if(length(TSSs) == 0) stop("No TSSs for FindMarkers.")
  object <- object[row.names(object) %in% TSSs, ]

  mat_p <- psi(object)
  mat_n <- counts(object)
  mat_n <- reshape2::melt(mat_n, value.name = "n")
  mat_p <- reshape2::melt(mat_p, value.name = "PSI")

  if(identical(mat_n[, 1:2], mat_p[, 1:2])) {
    mat <- cbind(mat_n, PSI = mat_p[, 3])
  } else {
    mat <- merge(mat_n, mat_p, by = c("Var1", "Var2"))
  }

  mat$groupBy <- plyr::mapvalues(as.character(mat$Var2), row.names(SummarizedExperiment::colData(object)), SummarizedExperiment::colData(object)[, groupBy])
  mat$gene_id <- plyr::mapvalues(as.character(mat$Var1), row.names(SummarizedExperiment::rowData(object)), SummarizedExperiment::rowData(object)[, "gene_id"])
  mat <- data.table::as.data.table(mat)
  data.table::setkey(mat, gene_id, Var2)

  N <- suppressWarnings(mat[, .(N = max(n/PSI, na.rm = T)), .(Var2, gene_id)])
  N[is.infinite(N), N := 0]
  mat <- data.table::merge.data.table(mat, N, by = c("gene_id", "Var2"))
  data.table::setnames(mat, "groupBy", groupBy)
  if(majorOnly) {
    mT <- data.table::as.data.table(as.data.frame(SummarizedExperiment::rowData(object)), keep.rownames = "ID")[, .SD[which.max(pro), ], gene_id][, ID]
    mat <- mat[Var1 %in% mT]
  }
  mat <- as.data.frame(mat)
  mat <- with(mat, split(mat, as.character(Var1)))
  res <- parallel::mclapply(mat, FUN = function(x) DTSS(x, groupBy, groups, group1, ratio1, ratio2, method, gs, deltaPSI), mc.cores = cores)
  data.table::as.data.table(do.call(rbind, do.call(c, res)))
}
