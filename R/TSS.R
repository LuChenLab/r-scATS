#' @name TSSCDF
#' @title TSS inference and quantification
#' @importFrom rtracklayer readGFFAsGRanges
#' @importFrom S4Vectors mcols
#' @importFrom pbapply pblapply
#' @importFrom parallel mclapply
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomeInfoDb keepSeqlevels
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
#' @param scDR Whether to estimate the degradation rate of single cell TSS.
#' @param MinscTSSReads Minimum reads of single cell TSS reads for single cell degradation rate estimation.
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

TSSCDF <- function (object, bam, sep = "_", genes, gtfFile, txdb = NULL, UTROnly = FALSE, MinGeneReads = 500, scDR = FALSE, MinscTSSReads = 5,
                    mapqFilter = 255, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = NA, Upstream = 1000,
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
  stopifnot(Upstream >= 0)

  if (!file.exists(gtfFile)) {
    stop("GTF file is not exists.")
  } else {
    if (is.null(txdb)) {
      txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(gtfFile)))
    } else {
      if (!is(txdb, "TxDb")) {
        txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(gtfFile)))
      }
    }
    seqs <- unique(do.call(c, lapply(bam, function(x) names(Rsamtools::scanBamHeader(x)[[1]]$target))))
    GenomicFeatures::isActiveSeq(txdb) <- GenomeInfoDb::seqlevels(txdb) %in% seqs
    gtf <- rtracklayer::readGFFAsGRanges(gtfFile)
    gtf <- GenomeInfoDb::keepSeqlevels(gtf, intersect(GenomeInfoDb::seqlevels(gtf), seqs), pruning.mode = "coarse")
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
  Tars <- GenomicRanges::promoters(Tars, upstream = Upstream, downstream = GenomicRanges::width(Tars))

  if (length(Tars) == 0) stop("There is no gene for ATS identification!")
  # save.image(file = "/mnt/raid61/Personal_data/tangchao/Temp/image.RData")
  if (verbose) {
    ATS_count <- pbapply::pblapply(X = seq_along(Tars),
                                   function(i) {
                                     tryCatch(gmmc(object = object, region = Tars[i],
                                                    bams = bam, gtfrange = gtf, UTR = fiveUTR, EBG = EBG, IBG = IBG,
                                                    UTROnly = UTROnly, MinGeneReads = MinGeneReads, scDR = scDR, MinscTSSReads = MinscTSSReads,
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
                                      tryCatch(gmmc(object = object, region = Tars[i],
                                                     bams = bam, gtfrange = gtf, UTR = fiveUTR, EBG = EBG, IBG = IBG,
                                                     UTROnly = UTROnly, MinGeneReads = MinGeneReads, scDR = scDR, MinscTSSReads = MinscTSSReads,
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
  if(scDR) {
    theta <- as.matrix(do.call(rbind, lapply(ATS_count, function(x) x$mat_theta)))
    alpha <- as.matrix(do.call(rbind, lapply(ATS_count, function(x) x$mat_alpha)))
  }
  rR <- as(with(paras, TSS), "GRanges")
  S4Vectors::mcols(rR) <- paras
  names(rR) <- with(paras, paste(gene_name, TSS, sep = "@"))
  if(scDR) {
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = count, psi = psi, theta = theta, alpha = alpha), rowRanges = rR, colData = slot(object, "meta.data"))
  } else {
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = count, psi = psi), rowRanges = rR, colData = slot(object, "meta.data"))
  }
  EBG <- EBG[names(EBG) %in% paras$gene_id]
  Tars <- subset(Tars, gene_id %in% paras$gene_id)
  parameters <- list(bam = bam, sep = sep, gtf = gtfFile, UTROnly = UTROnly, MinGeneReads = MinGeneReads, scDR = scDR, MinscTSSReads = MinscTSSReads, Upstream = Upstream,
                     mapqFilter = mapqFilter, isSecondaryAlignment = isSecondaryAlignment, isSupplementaryAlignment = isSupplementaryAlignment, isDuplicate = isDuplicate,
                     min.TSS.percent = min.TSS.percent, min.local.percent = min.local.percent, p.cutoff = p.cutoff, window = window, greedy = greedy, TargetRegions = Tars, EBG = EBG)
  scATSDataSet(se = se, reductions = slot(object, "reductions"), project.name = project.name, parameters = parameters)
}





#' @title FindMarkers
#' @name FindMarkers
#' @description Find markers of given group
#' @importFrom SummarizedExperiment colData
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment rowData
#' @importFrom data.table as.data.table
#' @importFrom parallel mclapply
#' @importFrom rlang %||%
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

#' @title FindMarkers
#' @name FindMarkers
#' @export

FindMarkersByTheta <- function(object, gene = NULL, cells = NULL, groupBy, group1 = NULL, group2 = NULL, cores = 1) {
  stopifnot(is(object, "scATSDataSet"))
  cells <- cells %||% colnames(object)
  stopifnot(all(is.element(cells, colnames(object))))
  gene <- gene %||% SummarizedExperiment::rowData(object)$gene_name
  gene <- intersect(gene, SummarizedExperiment::rowData(object)$gene_name)
  if(length(gene) == 0) return(NULL)
  stopifnot(all(gene %in% SummarizedExperiment::rowData(object)$gene_name))

  stopifnot(!is.null(groupBy))
  if(!is.null(groupBy)) {
    stopifnot(length(groupBy) == 1)
    stopifnot(groupBy %in% colnames(SummarizedExperiment::colData(object)))
  }

  if(is.null(group1) & !is.null(group2)) stop()
  if(!is.null(group1) & !is.null(group2)) {
    stopifnot(all(c(group1, group2) %in% SummarizedExperiment::colData(object)[, groupBy]))
    object <- object[, SummarizedExperiment::colData(object)[, groupBy] %in% c(group1, group2)]
  }

  grouplevels <- as.character(unique(SummarizedExperiment::colData(object)[, groupBy]))
  stopifnot(length(grouplevels) > 1)

  if(is.null(group1)) {
    groups <- lapply(grouplevels, function(x) {
      group1 = x; group2 = setdiff(grouplevels, x)
      grouplevelsFrom <- c(group1, group2)
      if(length(group1) == 1 & length(group2) == 1) {
        grouplevelsTo <- c(group1, group2)
      }
      if(length(group1) != 1 & length(group2) == 1) {
        grouplevelsTo <- c(rep("group1", length(group1)), group2)
      }
      if(length(group1) == 1 & length(group2) != 1) {
        grouplevelsTo <- c(group1, rep("group2", length(group2)))
      }
      if(length(group1) != 1 & length(group2) != 1) {
        grouplevelsTo <- c(rep("group1", length(group1)), rep("group2", length(group2)))
      }
      list(group1 = group1, group2 = group2, From = grouplevelsFrom, To = grouplevelsTo)
    })
  } else {
    grouplevelsFrom <- c(group1, group2)
    if(length(group1) == 1 & length(group2) == 1) {
      grouplevelsTo <- c(group1, group2)
    }
    if(length(group1) != 1 & length(group2) == 1) {
      grouplevelsTo <- c(rep("group1", length(group1)), group2)
    }
    if(length(group1) == 1 & length(group2) != 1) {
      grouplevelsTo <- c(group1, rep("group2", length(group2)))
    }
    if(length(group1) != 1 & length(group2) != 1) {
      grouplevelsTo <- c(rep("group1", length(group1)), rep("group2", length(group2)))
    }
    groups <- list(list(group1 = group1, group2 = group2, From = grouplevelsFrom, To = grouplevelsTo))
  }

  parallel::mclapply(gene, function(x) {
    para <- subset(SummarizedExperiment::rowData(object), gene_name == x)
    lapply(row.names(para), function(y) {
      mat_t <- data.table::data.table(CB = colnames(object), Theta = as.numeric(theta(object, cells, y)), groupBy = SummarizedExperiment::colData(object)[, groupBy])
      lapply(groups, function(z) {
        percent1 <- mat_t[groupBy %in% z$group1, mean(!is.na(Theta)) * 100]
        percent2 <- mat_t[groupBy %in% z$group2, mean(!is.na(Theta)) * 100]
        cell1 <- mat_t[groupBy %in% z$group1, sum(!is.na(Theta))]
        cell2 <- mat_t[groupBy %in% z$group2, sum(!is.na(Theta))]
        mat_t[, groupBy2 := plyr::mapvalues(groupBy, z$From, z$To, warn_missing = F)]
        p <- tryCatch(mat_t[, suppressWarnings(wilcox.test(Theta ~ groupBy2))$p.value], error = function(e) NA)
        data.table::data.table(group1 = paste(z$group1, collapse = ","), cell1, cell2, percent1, percent2, p)
      }) -> res
      data.table::data.table(TSS = y, do.call(rbind, res))
    }) -> res
    res <- data.table::data.table(gene = x, do.call(rbind, res))

    if(nrow(para) == 1) {
      res2 <- res[, .(gene, TSS, group1, theta1 = as.numeric(cell1 > 0), theta2 = as.numeric(cell2 > 0))]
    } else {
      res2 <- lapply(groups, function(z) {
        res2 <- ThetaByGroup(object = object, gene = x, groupBy = groupBy, group1 = z$group1, group2 = z$group2, cells = cells)
        if(nrow(res2) == 0) {
          res2 <- res[, .(gene, TSS, group1, theta1 = 0, theta2 = 0)]
        } else {
          res2 <- merge(res2[Group == "group1" | is.element(Group, z$group1), .(group1 = Group, TSS, theta1 = theta)],
                        res2[Group == "group2" | is.element(Group, z$group2), .(TSS, theta2 = theta)], by = "TSS", all = TRUE)
          res2 <- data.table::data.table(gene = x, res2)
          res2[, TSS := paste0(gene, "@", TSS)]
        }
      })
      res2 <- do.call(rbind, res2)
    }
    merge(res2, res, by = c("gene", "TSS", "group1"), all = TRUE)
  }, mc.cores = cores) -> res
  do.call(rbind, res)
}


#' @title ThetaByGroup
#' @name ThetaByGroup
#'
#' @param object An scATSDataSet from TSS.
#' @param gene the name of gene.
#' @param cells Which cells to use.
#' @param groupBy The column name of group in scATSDataSet object meta.data.
#' @param group1,group2 Only for given group.
#'
#' @export

ThetaByGroup <- function(object, gene, cells = NULL, groupBy, group1 = NULL, group2 = NULL) {
  cells <- cells %||% colnames(object)
  object <- object[, cells]
  stopifnot(length(gene) == 1)
  stopifnot(gene %in% SummarizedExperiment::rowData(object)$gene_name)
  stopifnot(!is.null(groupBy))
  if(!is.null(groupBy)) {
    stopifnot(length(groupBy) == 1)
    stopifnot(groupBy %in% colnames(SummarizedExperiment::colData(object)))
  }
  para <- subset(SummarizedExperiment::rowData(object), gene_name == gene)
  if(nrow(para) == 1) return(NULL)
  if(is.null(group1) & !is.null(group2)) stop()
  if(!is.null(group1)) {
    grouplevels <- as.character(unique(SummarizedExperiment::colData(object)[, groupBy]))
    stopifnot(all(is.element(group1, grouplevels)))
    if(is.null(group2)) {
      group2 <- setdiff(grouplevels, group1)
    } else {
      stopifnot(all(is.element(group1, grouplevels)))
      stopifnot(!any(group1 %in% group2))
    }
  }
  if(!is.null(group1)) {
    grouplevelsFrom <- c(group1, group2)
    if(length(group1) == 1 & length(group2) == 1) {
      grouplevelsTo <- c(group1, group2)
    }
    if(length(group1) != 1 & length(group2) == 1) {
      grouplevelsTo <- c(rep("group1", length(group1)), group2)
    }
    if(length(group1) == 1 & length(group2) != 1) {
      grouplevelsTo <- c(group1, rep("group2", length(group2)))
    }
    if(length(group1) != 1 & length(group2) != 1) {
      grouplevelsTo <- c(rep("group1", length(group1)), rep("group2", length(group2)))
    }
  }

  gene_i <- para$gene_id[1]
  parameters <- S4Vectors::metadata(object)[[2]]

  map0 <- loadbamgene(bam = parameters$bam, sep = parameters$sep, exonByGene = parameters$EBG,
                      gene = gene_i, region = subset(parameters$TargetRegions, gene_id == gene_i),
                      txdb = NULL, exonIntersect = TRUE, junctionIntersect = FALSE,
                      isSupplementaryAlignment = parameters$isSupplementaryAlignment,
                      isSecondaryAlignment = parameters$isSecondaryAlignment,
                      isDuplicate = parameters$isDuplicate, mapqFilter = parameters$mapqFilter)
  map0 <- map0[S4Vectors::mcols(map0)$CB %in% cells]
  si <- as.character(BiocGenerics::strand(subset(parameters$TargetRegions, gene_id == gene_i)))

  if(si == "+") {
    Pos <- data.table::data.table(pos = IRanges::start(IRanges::ranges(map0)), CB = S4Vectors::mcols(map0)$CB)
  } else {
    Pos <- data.table::data.table(pos = IRanges::end(IRanges::ranges(map0)), CB = S4Vectors::mcols(map0)$CB)
  }

  Pos <- merge(Pos, data.table::data.table(CB = row.names(SummarizedExperiment::colData(object)), groupBy = SummarizedExperiment::colData(object)[, groupBy]), by = "CB")
  if(!is.null(group1)) {
    Pos <- Pos[groupBy %in% grouplevelsFrom]
    Pos[, groupBy := plyr::mapvalues(groupBy, grouplevelsFrom, grouplevelsTo, warn_missing = FALSE)]
  }

  TsPos <- lapply(seq_len(nrow(para)), function(x) {
    data.table::data.table(TSS = para$TSS[x], Pos[pos %inrange% as.data.frame(as(para$Region[x], "IRanges"))])
  })
  TsPos <- na.omit(do.call(rbind, TsPos))

  if(si == "+") {
    TsPos[, pos2 := pos]
    for (i in seq_len(nrow(para))) {
      bas <- IRanges::start(as(para[i, "Region"], "IRanges")):IRanges::end(as(para[i, "Region"], "IRanges"))
      TsPos[, pos2 := plyr::mapvalues(pos2, bas, pmax(bas, GenomicRanges::start(as(para[i, "TSS"], "GRanges"))), warn_missing = FALSE)] # Be careful for +/-
    }
    TsPos[, TSS := factor(TSS, levels = TsPos[, min(pos2), TSS][order(V1, decreasing = F), TSS])] # Be careful for +/-
  } else {
    TsPos[, pos2 := pos]
    for (i in seq_len(nrow(para))) {
      bas <- IRanges::end(as(para[i, "Region"], "IRanges")):IRanges::start(as(para[i, "Region"], "IRanges"))
      TsPos[, pos2 := plyr::mapvalues(pos2, bas, pmin(bas, GenomicRanges::start(as(para[i, "TSS"], "GRanges"))), warn_missing = FALSE)] # Be careful for +/-
    }
    TsPos[, TSS := factor(TSS, levels = TsPos[, max(pos2), TSS][order(V1, decreasing = T), TSS])] # Be careful for +/-
  }

  if(si == "+") {
    DRT <- TsPos[, .(count = .N), .(TSS, start = pos2, index = as.numeric(TSS), groupBy)] # Be careful for +/-
  } else {
    DRT <- TsPos[, .(count = .N), .(TSS, start = -pos2, index = as.numeric(TSS), groupBy)] # Be careful for +/-
  }
  DRT$start <- DRT$start - min(DRT$start) + 1
  DRT[, end := start]
  DRT <- DRT[order(groupBy, index, start)]
  DRT <- split(DRT, DRT$groupBy)
  DRT <- lapply(DRT, function(x) {x$count <- cumsum(x$count); return(x)})
  DRR <- lapply(DRT, function(x) tryCatch(DR(x = x, minreads = parameters$MinscTSSReads), error = function(e) NULL))
  DRR <- DRR[!mapply(is.null, DRR)]
  DRR <- data.table::data.table(Group = rep(names(DRR), mapply(nrow, DRR)), do.call(rbind, DRR))
  return(DRR)
}

