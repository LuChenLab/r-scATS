#' @title InternalFunctions
#' @description Internal function for FindMarkers
#' @importFrom IRanges IRanges
#' @importFrom IRanges start
#' @importFrom IRanges end
#'


#' @title TSS count from GMM
#' @description TSS count from GMM
#'
#' @importFrom rtracklayer readGFFAsGRanges
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicAlignments findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols
#' @importFrom GenomicAlignments first
#' @importFrom GenomicAlignments second
#' @importFrom BiocGenerics strand
#' @importFrom data.table data.table setnames
#' @importFrom IRanges start
#' @importFrom IRanges end
#' @importFrom IRanges ranges
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges start
#' @importFrom plyr mapvalues
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table dcast
#' @importFrom Seurat FetchData
#' @importFrom GenomeInfoDb seqnames
#' @importFrom mclust densityMclust
#' @importFrom mclust priorControl
#' @importFrom mixtools normalmixEM
#' @importFrom R.utils withTimeout
#' @importFrom GenomicRanges promoters
#' @importFrom GenomicRanges reduce
#' @importFrom utils capture.output
#'

#' @title DTSS
#' @description Differentially expressed TSS
#' @importFrom plyr mapvalues
#' @importFrom SummarizedExperiment colData
#' @importFrom VGAM vglm
#' @importFrom lmtest waldtest
DTSS <- function(x, groupBy, groups, group1, ratio1, ratio2, method, gs, deltaPSI = 0) {
  lapply(groups, function(g) {
    x[groupBy] <- plyr::mapvalues(x[, groupBy], gs, g)
    # mati <- data.table::as.data.table(subset(x, !is.na(PSI)))
    mati <- data.table::as.data.table(x)
    mati[is.na(n), n := 0]

    if(length(unique(mati[, groupBy, with = FALSE][[1]])) < 2) {
      return(NULL)
    } else {
      colnames(mati)[colnames(mati) == groupBy] <- "groupBy"
      if(is.null(group1)) {
        G1 <- g[g != "G0"]
        G2 <- "Other"
        n1 <- sum(with(mati, groupBy != "G0" & n > 0))
        n2 <- sum(with(mati, groupBy == "G0" & n > 0))
        N1 <- sum(with(mati, groupBy != "G0" & N > 0))
        N2 <- sum(with(mati, groupBy == "G0" & N > 0))
        Cells1 <- mati[groupBy != "G0", .N]
        Cells2 <- mati[groupBy == "G0", .N]
      } else {
        G1 <- group1
        G2 <- setdiff(g, G1)
        G2 <- G2[!is.na(G2)]
        n1 <- sum(with(mati, groupBy == G1 & n > 0))
        n2 <- sum(with(mati, groupBy != G1 & n > 0))
        N1 <- sum(with(mati, groupBy == G1 & N > 0))
        N2 <- sum(with(mati, groupBy != G1 & N > 0))
        Cells1 <- mati[groupBy == G1, .N]
        Cells2 <- mati[groupBy != G1, .N]
      }
      if(is.na(n1/N1) & is.na(n2/N2)) return(NULL)
      if(ifelse(is.na(n1/N1), 1, n1/N1) < ratio1 | ifelse(is.na(n2/N2), 1, n2/N2) < ratio2) {
        return(NULL)
      } else {
        PSI1 = mati[!is.na(PSI) & groupBy == G1, mean(PSI)]
        PSI2 = mati[!is.na(PSI) & groupBy != G1, mean(PSI)]
        PseudobulkPSI1 = mati[!is.na(PSI) & groupBy == G1, sum(n) / sum(N)]
        PseudobulkPSI2 = mati[!is.na(PSI) & groupBy != G1, sum(n) / sum(N)]

        if(ifelse(is.na(abs(PSI1 - PSI2)), max(c(PSI1, PSI2), na.rm = T), abs(PSI1 - PSI2)) < deltaPSI & ifelse(is.na(abs(PseudobulkPSI1 - PseudobulkPSI2)), max(c(PseudobulkPSI1, PseudobulkPSI2), na.rm = T), abs(PseudobulkPSI1 - PseudobulkPSI2)) < deltaPSI) {
          return(NULL)
        } else {
          if(mati[!is.na(PSI), sd(PSI)] == 0 | is.na(n1/N1) | is.na(n2/N2)) {
            if(is.na(n1/N1) | is.na(n2/N2)) {
              if(is.null(method)) {
                return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                  PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                  wald.test = 0, wilcox.test = 0, prop.test = 0))
              } else {
                if(method == "wilcox.test") {
                  return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                    PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                    wilcox.test = 0))
                }

                if(method == "prop.test") {
                  return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                    PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                    prop.test = 0))
                }

                if(method == "wald.test") {
                  return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                    PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                    wald.test = 0))
                }
              }
            } else {
              if(is.null(method)) {
                return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                  PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                  wald.test = NA, wilcox.test = NA, prop.test = NA))
              } else {
                if(method == "wilcox.test") {
                  return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                    PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                    wilcox.test = NA))
                }

                if(method == "prop.test") {
                  return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                    PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                    prop.test = NA))
                }

                if(method == "wald.test") {
                  return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                    PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                    wald.test = NA))
                }
              }
            }
          } else {
            if(is.null(method)) {
              prop.test <- prop.test(matrix(c(mati[groupBy == G1, sum(n)],
                                              mati[groupBy == G1, sum(N) - sum(n)],
                                              mati[groupBy != G1, sum(n)],
                                              mati[groupBy != G1, sum(N) - sum(n)]), ncol = 2, byrow = T))
              wilcox.test <- mati[!is.na(PSI), suppressWarnings(wilcox.test(PSI ~ groupBy))]
              m0 <- tryCatch(suppressWarnings(VGAM::vglm(cbind(n, N - n) ~ 1, "betabinomial", data = mati[!is.na(PSI)], trace = FALSE, maxit = 100)), error = function(e) NA)
              m1 <- tryCatch(suppressWarnings(VGAM::vglm(cbind(n, N - n) ~ groupBy, "betabinomial", data = mati[!is.na(PSI)], trace = FALSE, maxit = 100)), error = function(e) NA)
              wald.test <- tryCatch(suppressWarnings(lmtest::waldtest(m0, m1))$`Pr(>Chisq)`[2], error = function(e) NA)
              return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                wald.test, wilcox.test = wilcox.test$p.value, prop.test = prop.test$p.value))
            } else {
              if(method == "wilcox.test") {
                wilcox.test <- mati[!is.na(PSI), suppressWarnings(wilcox.test(PSI ~ groupBy))]
                return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                  PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                  wilcox.test = wilcox.test$p.value))
              }
              if(method == "prop.test") {
                prop.test <- prop.test(matrix(c(mati[groupBy == G1, sum(n)],
                                                mati[groupBy == G1, sum(N) - sum(n)],
                                                mati[groupBy != G1, sum(n)],
                                                mati[groupBy != G1, sum(N) - sum(n)]), ncol = 2, byrow = T))
                return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                  PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                  prop.test = prop.test$p.value))
              }
              if(method == "wald.test") {
                m0 <- tryCatch(suppressWarnings(VGAM::vglm(cbind(n, N - n) ~ 1, "betabinomial", data = mati[!is.na(PSI)], trace = FALSE, maxit = 100)), error = function(e) NA)
                m1 <- tryCatch(suppressWarnings(VGAM::vglm(cbind(n, N - n) ~ groupBy, "betabinomial", data = mati[!is.na(PSI)], trace = FALSE, maxit = 100)), error = function(e) NA)
                wald.test <- tryCatch(suppressWarnings(lmtest::waldtest(m0, m1))$`Pr(>Chisq)`[2], error = function(e) NA)
                return(data.frame(TSS = unique(as.character(mati$Var1)), G1, G2, n1, n2, N1, N2, Cells1, Cells2,
                                  PSI1 = PSI1, PSI2 = PSI2, PseudobulkPSI1 = PseudobulkPSI1, PseudobulkPSI2 = PseudobulkPSI2,
                                  wald.test))
              }
            }
          }
        }
      }
    }
  })
}

#' @title InternalFunctions
#' @description Internal function for TSSPlot, read bam reads and calculate start site and coverage for
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicAlignments findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols
#' @importFrom GenomicAlignments first
#' @importFrom GenomicAlignments second
#' @importFrom BiocGenerics strand
#' @importFrom GenomicAlignments coverage
#' @importFrom S4Vectors runValue
#'
rss <- function(region, bam, EBG, ExonicReadsOnly = FALSE, mapqFilter = 255,
                isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = NA, sep = "_") {
  param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = isSupplementaryAlignment,
                                                                 isPaired = TRUE,
                                                                 isDuplicate = isDuplicate,
                                                                 isNotPassingQualityControls = FALSE,
                                                                 isSecondaryAlignment = isSecondaryAlignment),
                                   tag = "CB", mapqFilter = mapqFilter, which = region)
  if(length(bam) == 1) {
    bami <- paste0(tempfile(), ".bam")
    sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, bam, gsub(":[-\\+]", "", as.character(region)))
    system(command = sh); system(command = paste("samtools index", bami))
    map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bami, param = param, use.names = TRUE, strandMode = 1))
    R2 <- GenomicAlignments::second(map0)
    GenomicRanges::strand(R2) <- GenomicRanges::strand(map0)
    map0 <- map0[S4Vectors::queryHits(GenomicAlignments::findOverlaps(R2, region, type = "within", ignore.strand = FALSE))]
    map0 <- map0[S4Vectors::mcols(GenomicAlignments::first(map0))$CB == S4Vectors::mcols(GenomicAlignments::second(map0))$CB & !is.na(S4Vectors::mcols(GenomicAlignments::first(map0))$CB)]
    if(ExonicReadsOnly) {
      R2 <- GenomicAlignments::second(map0)
      GenomicRanges::strand(R2) <- GenomicRanges::strand(map0)
      R2_M <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = GenomicAlignments::cigar(R2), ops = "M", pos = GenomicAlignments::start(R2))
      names(R2_M) <- names(R2)
      R2Us <- unique(names(unlist(R2_M))[S4Vectors::queryHits(IRanges::findOverlaps(unlist(R2_M), IRanges::ranges(EBG[[region$gene_id]]), type = "within"))])
      map0 <- map0[names(map0) %in% R2Us]
    }
    if(length(map0) > 0) {
      S4Vectors::mcols(map0)$CB <- S4Vectors::mcols(GenomicAlignments::first(map0))$CB
    }
    system(command = paste0("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
  } else {
    map0 <- lapply(bam, function(x) {
      bami <- paste0(tempfile(), "_", basename(x))
      sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, x, gsub(":[-\\+]", "", as.character(region)))
      system(command = sh); system(command = paste("samtools index", bami))
      map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bami, param = param, use.names = TRUE, strandMode = 1))
      # map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(x, param = param, use.names = TRUE, strandMode = 1))
      R2 <- GenomicAlignments::second(map0)
      GenomicRanges::strand(R2) <- GenomicRanges::strand(map0)
      map0 <- map0[S4Vectors::queryHits(GenomicAlignments::findOverlaps(R2, region, type = "within", ignore.strand = FALSE))]
      map0 <- map0[S4Vectors::mcols(GenomicAlignments::first(map0))$CB == S4Vectors::mcols(GenomicAlignments::second(map0))$CB & !is.na(S4Vectors::mcols(GenomicAlignments::first(map0))$CB)]
      if(ExonicReadsOnly) {
        R2 <- GenomicAlignments::second(map0)
        GenomicRanges::strand(R2) <- GenomicRanges::strand(map0)
        R2_M <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = GenomicAlignments::cigar(R2), ops = "M", pos = GenomicAlignments::start(R2))
        names(R2_M) <- names(R2)
        R2Us <- unique(names(unlist(R2_M))[S4Vectors::queryHits(IRanges::findOverlaps(unlist(R2_M), IRanges::ranges(EBG[[region$gene_id]]), type = "within"))])
        map0 <- map0[names(map0) %in% R2Us]
      }
      if(length(map0) > 0) {
        S4Vectors::mcols(map0)$CB <- paste(S4Vectors::mcols(GenomicAlignments::first(map0))$CB, gsub(".bam", "", basename(x)), sep = sep)
      }
      system(command = paste("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
      return(map0)
    })
    map0 <- do.call(c, map0)
  }
  if(length(map0) == 0) return(NULL)

  if(as.character(BiocGenerics::strand(region)) == "+") {
    Pos <- data.table::data.table(pos = IRanges::start(IRanges::ranges(map0)), CB = S4Vectors::mcols(map0)$CB)
    cov <- GenomicAlignments::coverage(GenomicAlignments::first(map0))[[as.character(S4Vectors::runValue(GenomicAlignments::seqnames(map0)))]]
  } else {
    Pos <- data.table::data.table(pos = IRanges::end(IRanges::ranges(map0)), CB = S4Vectors::mcols(map0)$CB)
    cov <- GenomicAlignments::coverage(GenomicAlignments::first(map0))[[as.character(S4Vectors::runValue(GenomicAlignments::seqnames(map0)))]]
  }
  cov <- data.table::data.table(pos = Pos[, min(pos)]:Pos[, max(pos)],
                                D = as.numeric(cov[Pos[, min(pos)]:Pos[, max(pos)]]))
  return(list(cov, Pos))
}




#' @title InternalFunctions
#' @description Internal function for TSSCDF
#' @importFrom IRanges IRanges countOverlaps
#' @importFrom S4Vectors subjectHits mcols
gmmc <- function(object, region, bams, gtfrange, UTR, EBG, IBG, UTROnly = FALSE, greedy = TRUE,
                  MinGeneReads = 1000, min.TSS.percent = 5, min.local.percent = 1, p.cutoff = 0.01, window = 10, scDR = FALSE, MinscTSSReads = 5,
                  mapqFilter = 255, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = FALSE, sep = "_") {
  CellsTU <- as.character(slot(object, "Cell"))
  map0 <- loadbamgene(bam = bams, gene = with(as.data.frame(region), gene_id), region = region,
                      txdb = NULL, exonByGene = EBG,
                      exonIntersect = TRUE, junctionIntersect = FALSE,
                      isSupplementaryAlignment = isSupplementaryAlignment,
                      isSecondaryAlignment = isSecondaryAlignment,
                      isDuplicate = isDuplicate, mapqFilter = mapqFilter, sep = sep)
  if(is.null(map0)) return(NULL)
  map0 <- map0[S4Vectors::mcols(map0)$CB %in% CellsTU]
  if(length(map0) < MinGeneReads) {
    warning("Reads < MinGeneReads.")
    return(NULL)
  }

  if(as.character(BiocGenerics::strand(region)) == "+") {
    Pos <- data.table::data.table(pos = IRanges::start(IRanges::ranges(GenomicAlignments::first(map0))), CB = S4Vectors::mcols(map0)$CB)
  } else {
    Pos <- data.table::data.table(pos = IRanges::end(IRanges::ranges(GenomicAlignments::first(map0))), CB = S4Vectors::mcols(map0)$CB)
  }
  rm("map0")
  if(nrow(Pos) < MinGeneReads)  {
    warning("Reads < MinGeneReads")
    return(NULL)
  }

  if(is.null(UTROnly)) {
    UTROnly <- Pos[pos %in% Pos[, .N, pos][N > 1, pos], mean(data.table::inrange(pos, start(EBG[[with(as.data.frame(region), gene_id)]]), end(EBG[[with(as.data.frame(region), gene_id)]])))] >= 0.8
  }
  if(UTROnly) {
    utr <- UTR[UTR$gene_name == region$gene_name]
    utr5 <- GenomicRanges::reduce(utr)
    Pos <- Pos[S4Vectors::queryHits(IRanges::findOverlaps(Pos[, IRanges::IRanges(pos)], IRanges::ranges(utr5), maxgap = window)), ]
  }

  if(nrow(Pos) < MinGeneReads)  {
    warning("Reads < MinGeneReads")
    return(NULL)
  }
  tsss <- unique(GenomicRanges::start(GenomicRanges::promoters(gtfrange[gtfrange$type == "transcript" & gtfrange$gene_name == region$gene_name], upstream = 0, downstream = 1)))
  if(Pos[, sd(pos)] == 0) {
    paras <- data.table::data.table(gene_id = region$gene_id,
                                    gene_name = region$gene_name,
                                    TSS = paste(as.character(GenomeInfoDb::seqnames(region)), Pos[, unique(pos)], as.character(GenomicRanges::strand(region)), sep = ":"),
                                    Region = paste0(Pos[, unique(pos)], "-", Pos[, unique(pos)]), PSI = 1, Percent = 1,
                                    Annotated = ifelse(any(is.element(tsss, Pos[, unique(pos)])), Pos[, unique(pos)], NA),
                                    Greedy = FALSE, AllReads = Pos[ , .N], beta = 1, alpha = 0, theta = 1)
    Pos$TSS <- paras$TSS
    mat <- Pos[, .(count = .N), .(CB, TSS)]
    data.table::setkey(mat, CB)
    if(any(duplicated(mat[CellsTU, ][, .(TSS, CB)]))) return(NULL)
    mat_count <- data.table::dcast(mat[CellsTU, ], TSS ~ CB, value.var = "count")
    mat_count <- mat_count[!is.na(TSS)]
    data.table::setkey(mat_count, TSS)

    if(!identical(mat_count$TSS, paras$TSS)) mat_count <- mat_count[paras$TSS, ]
    mat_count[, TSS := NULL]

    if(!identical(colnames(mat_count), CellsTU)) mat_count <- mat_count[, CellsTU, with = F]

    if(paras[, .N] > 1) {
      mat_psi <- apply(mat_count, 2, myproportions)
    } else {
      mat_psi <- matrix(apply(mat_count, 2, myproportions), nrow = 1)
      colnames(mat_psi) <- colnames(mat_count)
    }
    if(scDR) {
      mat_alpha <- matrix(data = 1, nrow = nrow(mat_psi), ncol = ncol(mat_psi), dimnames = list(rownames(mat_psi), colnames(mat_psi)))
      res <- list(mat_count = mat_count, mat_psi = mat_psi, mat_theta = mat_psi, mat_alpha = mat_alpha, paras = paras)
    } else {
      res <- list(mat_count = mat_count, mat_psi = mat_psi, paras = paras)
    }
  } else {
    if(as.character(GenomicRanges::strand(region)) == "+") {
      Ts <- CliffSite(x = Pos$pos, direction = "start", introns = IBG[names(IBG) == with(as.data.frame(region), gene_id)], RefTSS = tsss, plots = 3,
                      minreads = MinGeneReads, min.TSS.percent = min.TSS.percent, min.local.percent = min.local.percent, p.cutoff = p.cutoff, window = window, greedy = greedy)
    } else {
      Ts <- CliffSite(x = Pos$pos, direction = "end", introns = IBG[names(IBG) == with(as.data.frame(region), gene_id)], RefTSS = tsss, plots = 3,
                      minreads = MinGeneReads, min.TSS.percent = min.TSS.percent, min.local.percent = min.local.percent, p.cutoff = p.cutoff, window = window, greedy = greedy)
    }
    if(length(Ts) == 0) return(NULL)

    TsPos <- lapply(seq_along(Ts), function(x) {
      Pos[pos %in% with(as.data.frame(Ts[x]), start:end)]
    })
    Ts <- as.data.table(Ts)
    paras <- data.table::data.table(gene_id = region$gene_id,
                                    gene_name = region$gene_name,
                                    TSS = paste(as.character(GenomeInfoDb::seqnames(region)), with(Ts, SS), as.character(BiocGenerics::strand(region)), sep = ":"),
                                    Region = with(Ts, paste(start, end, sep = "-")),
                                    PSI = with(Ts, PSI), Percent = with(Ts, Percent),
                                    Annotated = with(Ts, Annotated), Greedy = with(Ts, Greedy),
                                    AllReads = Pos[ , .N], beta = with(Ts, beta))
    TsPos <- data.table::data.table(TSS = rep(paras[, TSS], mapply(nrow, TsPos)), do.call(rbind, TsPos))

    if(as.character(GenomicRanges::strand(region)) == "+") {
      TsPos[, pos2 := pos]
      for (i in seq_len(nrow(Ts))) {
        TsPos[, pos2 := plyr::mapvalues(pos2, Ts[i, start]:Ts[i, end], pmax(Ts[i, start]:Ts[i, end], Ts[i, SS]), warn_missing = FALSE)] # Be careful for +/-
      }
      TsPos[, TSS := factor(TSS, levels = TsPos[, min(pos2), TSS][order(V1, decreasing = F), TSS])] # Be careful for +/-
      DRT <- TsPos[, .(count = .N), .(TSS, start = pos2, index = as.numeric(TSS))] # Be careful for +/-
    } else {
      TsPos[, pos2 := pos]
      for (i in seq_len(nrow(Ts))) {
        TsPos[, pos2 := plyr::mapvalues(pos2, Ts[i, end]:Ts[i, start], pmin(Ts[i, end]:Ts[i, start], Ts[i, SS]), warn_missing = FALSE)] # Be careful for +/-
      }
      TsPos[, TSS := factor(TSS, levels = TsPos[, max(pos2), TSS][order(V1, decreasing = T), TSS])] # Be careful for +/-
      DRT <- TsPos[, .(count = .N), .(TSS, start = -pos2, index = as.numeric(TSS))] # Be careful for +/-
    }

    DRT$start <- DRT$start - min(DRT$start) + 1
    DRT[, end := start]
    DRT <- DRT[order(index, start)]
    DRT$count <- cumsum(DRT$count)
    DRR <- DR(x = DRT)
    if(!is.null(DRR)) {
      paras <- merge(paras, DRR, by = "TSS", all.x = TRUE)
    } else {
      paras <- data.table(paras, alpha = NA, theta = NA)
    }
    data.table::setkey(paras, TSS)

    if(scDR & paras[, .N] >= 1) {
      if(as.character(GenomicRanges::strand(region)) == "+") {
        DRT <- TsPos[, .(count = .N), .(TSS, start = pos2, index = as.numeric(TSS), CB)] # Be careful for +/-
      } else {
        DRT <- TsPos[, .(count = .N), .(TSS, start = -pos2, index = as.numeric(TSS), CB)] # Be careful for +/-
      }
      DRT$start <- DRT$start - min(DRT$start) + 1
      DRT[, end := start]
      DRT <- DRT[order(CB, index, start)]
      DRT <- split(DRT, DRT$CB)
      DRT <- lapply(DRT, function(x) {x$count <- cumsum(x$count); return(x)})
      DRR <- lapply(DRT, function(x) DR(x = x, minreads = MinscTSSReads))
      DRR <- DRR[!mapply(is.null, DRR)]
      DRR <- data.table::data.table(CB = rep(names(DRR), mapply(nrow, DRR)), do.call(rbind, DRR))

      mat_alpha <- data.table::dcast(DRR, TSS ~ CB, value.var = "alpha")
      mat_alpha <- mat_alpha[!is.na(TSS)]
      data.table::setkey(mat_alpha, TSS)
      if(!identical(mat_alpha$TSS, paras$TSS)) mat_alpha <- mat_alpha[paras$TSS, ]
      mat_alpha[, TSS := NULL]
      mat_alpha <- data.table::as.data.table(t(mat_alpha), keep.rownames = "CB")
      data.table::setkey(mat_alpha, CB)
      if(!identical(mat_alpha[, CB], CellsTU)) mat_alpha <- mat_alpha[CellsTU, ]
      mat_alpha <- as.data.table(t(data.frame(mat_alpha[, -1], row.names = mat_alpha[[1]])))

      mat_theta <- data.table::dcast(DRR, TSS ~ CB, value.var = "theta")
      mat_theta <- mat_theta[!is.na(TSS)]
      data.table::setkey(mat_theta, TSS)
      if(!identical(mat_theta$TSS, paras$TSS)) mat_theta <- mat_theta[paras$TSS, ]
      mat_theta[, TSS := NULL]
      mat_theta <- data.table::as.data.table(t(mat_theta), keep.rownames = "CB")
      data.table::setkey(mat_theta, CB)
      if(!identical(mat_theta[, CB], CellsTU)) mat_theta <- mat_theta[CellsTU, ]
      mat_theta <- as.data.table(t(data.frame(mat_theta[, -1], row.names = mat_theta[[1]])))
    }
    mat <- TsPos[, .(count = .N), .(CB, TSS)]
    data.table::setkey(mat, CB)
    if(any(duplicated(mat[CellsTU, ][, .(TSS, CB)]))) return(NULL)

    mat_count <- data.table::dcast(mat[CellsTU, ], TSS ~ CB, value.var = "count")
    mat_count <- mat_count[!is.na(TSS)]
    data.table::setkey(mat_count, TSS)

    if(!identical(mat_count$TSS, paras$TSS)) mat_count <- mat_count[paras$TSS, ]
    mat_count[, TSS := NULL]

    if(!identical(colnames(mat_count), CellsTU)) mat_count <- mat_count[, CellsTU, with = F]
    if(paras[, .N] > 1) {
      mat_psi <- apply(mat_count, 2, myproportions)
    } else {
      mat_psi <- matrix(apply(mat_count, 2, myproportions), nrow = 1)
      colnames(mat_psi) <- colnames(mat_count)
    }

    if(scDR) {
      if(paras[, .N] > 1) {
        res <- list(mat_count = mat_count, mat_psi = mat_psi, mat_theta = mat_theta, mat_alpha = mat_alpha, paras = paras)
      } else {
        res <- list(mat_count = mat_count, mat_psi = mat_psi, mat_theta = mat_psi, mat_alpha = mat_alpha, paras = paras)
      }
    } else {
      res <- list(mat_count = mat_count, mat_psi = mat_psi, paras = paras)
    }
  }
  return(res)
}


#' @title InternalFunctions
#' @description Internal function for gmmc
#' @importFrom plyr mapvalues
#' @importFrom MASS fitdistr
#' @importFrom data.table `%inrange%` inrange
#' @importFrom BiocGenerics sort

CliffSite <- function(x, introns = NULL, direction = "start", minreads = 100, min.TSS.percent = 5, min.local.percent = 1, p.cutoff = 0.01,
                      window = 6, RefTSS = NULL, greedy = FALSE, plots = 0) {
  if(length(x) < minreads) {
    message("The number of reads is smaller than minreads")
    return(NULL)
  }
  stopifnot(is.element(direction, c("start", "end")))
  if(!is.null(introns)) {
    sjsite <- if(direction == "start") end(introns) + 1 else start(introns) - 1
  } else {
    sjsite <- NULL
  }
  From <- sort(unique(c(x, RefTSS, sjsite)))
  From <- c((min(From) - window):(min(From) - 1), From, (max(From) + 1):(max(From) + window))
  To <- if(direction == "start") seq_along(From) else rev(seq_along(From))

  RefTSS <- plyr::mapvalues(RefTSS, From, To, warn_missing = FALSE)
  sjsite <- plyr::mapvalues(sjsite, From, To, warn_missing = FALSE)
  x <- plyr::mapvalues(x, From, To, warn_missing = FALSE)

  if(length(x) == 0) return(NULL)
  if(length(RefTSS) == 0) RefTSS <- NULL

  set.seed(length(x))
  F0 <- ecdf(sample(min(x):max(x), size = length(x), replace = TRUE))
  F1 <- ecdf(x)
  L1 <- data.table::data.table(pos = min(x):max(x), CDR0 = F0(min(x):max(x)), CDR1 = F1(min(x):max(x)))
  L1 <- merge(data.table::data.table(pos = x)[, .N, pos], L1, by = "pos", all.y = T)[order(pos)]
  data.table::setkey(L1, pos)
  L1[is.na(N), N := 0]

  if(L1[, .N] > 1) {
    L1$CDR11 <- smthgaussian(L1$CDR1, window = pmin(window, L1[, .N - 1]), tails = TRUE)
    L1[, D := CDR1 - CDR11]
    L1$d0 <- L1[, splinefun(pos, CDR0, "monoH.FC")(pos, deriv = 1)]
    L1$d1 <- L1[, splinefun(pos, CDR11, "monoH.FC")(pos, deriv = 1)]
    # L1$d1 <- smthgaussian(L1$d1, window = pmin(window, L1[, .N - 1]), tails = TRUE)

    # Peaks
    m0 <- InflexionPoint(x = L1$pos, y = L1$d1)
    m1 <- m0[P == "T"]; m0 <- m0[P == "B", ]
    m1$start <- mapply(m1$pos, FUN = function (m) m0[pos < m, ifelse(.N == 0, min(x), max(end))])
    m1$end <- mapply(m1$pos, FUN = function (m) m0[pos > m, ifelse(.N == 0, max(x), min(start))])
    m1[, P := (F1(end) - F1(start)) * 100]
    m1 <- m1[P > 0]

    # RefTSS in region
    if(length(RefTSS) > 0) RefTSS <- RefTSS[RefTSS %inrange% m1[, .(start, end)]]
    if(length(RefTSS) > 0) {
      if(any(RefTSS %inrange% m1[, .(start, end)])) {
        ol <- IRanges::findOverlaps(IRanges::IRanges(RefTSS, RefTSS), m1[, IRanges::IRanges(start, end)])
        ol <- with(as.data.frame(ol), data.table::data.table(RefTSS = RefTSS[queryHits], pos = m1[subjectHits, pos]))
        ol <- ol[, .(RefTSS = RefTSS[which.min(abs(RefTSS - pos))]), pos]
        m1 <- merge(m1, ol, by = "pos", all.x = TRUE)
      } else {
        m1[, RefTSS := NA]
      }
    } else {
      m1[, RefTSS := NA]
    }

    # Candidate 1
    mod0 <- MASS::fitdistr(L1[!is.na(d0), d0], "normal")
    ci <- qnorm(1 - p.cutoff, mod0$estimate[1], mod0$estimate[2])
    cd1 <- data.table::as.data.table(L1[d1 > ci, IRanges::IRanges(pos)]); cd1[, width := NULL]
    if(nrow(cd1) > 0) {
      cd1 <- m1[pos %inrange% cd1, .(start, end, Ref = !is.na(RefTSS))]
      cd1[, P := (F1(end) - F1(start)) * 100]
      cd1 <- cd1[Ref | P >= min.TSS.percent]
    }
    if(nrow(cd1) > 0) cd1 <- cd1[mapply(function(i) L1[pos %inrange% cd1[i], max(abs(D))] , seq_len(nrow(cd1))) >= min.local.percent / 100]
    if(nrow(cd1) > 0) cd1 <- cd1[mapply(seq_len(nrow(cd1)), FUN = function(i) m1[!is.na(RefTSS) | P >= min.TSS.percent, any(pos %inrange% cd1[i])])] # Has significant peak
    if(nrow(cd1) > 0) if(cd1[, sum(P)] <= 50 & cd1[, max(P)] <= 10) cd1 <- cd1[P == 100]
    cd1$SS <- mapply(function(i) {
      if(m1[pos %inrange% cd1[i], any(P >= min.TSS.percent)]) {
        L1[pos %inrange% cd1[i]][which.max(d1), pos]
      } else {
        NA
      }
    }, seq_len(nrow(cd1)))

    cd1$Ref <- mapply(function(i) {
      if(cd1[i, Ref]) {
        # m1[pos %inrange% cd1[i] & !is.na(RefTSS)][which.max(P), RefTSS]
        L1[pos %inrange% cd1[i] & pos %in% RefTSS][which.max(d1), pos]
      } else {
        NA
      }
    }, seq_len(nrow(cd1)))

    cd1 <- cd1[!is.na(SS)]
    if(length(x) < minreads & nrow(cd1) > 0) cd1 <- cd1[P > 1]

    # First TSS
    cd1[, Greedy := FALSE]
    if(greedy == TRUE & m1[, any(!is.na(RefTSS))]) {
      if(nrow(cd1) > 0) {
        if(any(m1[!is.na(RefTSS), end] < cd1[, min(start)])) {
          cd1 <- rbind(m1[!is.na(RefTSS) & end < cd1[, min(start)], .(start, end, P, Ref = RefTSS, SS = RefTSS, Greedy = TRUE)][which.min(SS)], cd1)
        }
      } else {
        cd1 <- m1[!is.na(RefTSS), .(start, end, P, Ref = RefTSS, SS = RefTSS, Greedy = TRUE)][which.min(SS)]
      }
    }

    # remove splice junction sites
    if(!is.null(sjsite)) {
      if(nrow(cd1) > 0) cd1 <- cd1[with(as.data.frame(IRanges::distanceToNearest(cd1[, IRanges::IRanges(SS, SS)], IRanges::IRanges(sjsite))), distance) >= 3]
    }
    if(nrow(cd1) == 0) {
      solution <- NULL
    } else {
      solution <- cd1[, IRanges::IRanges(start, end, SS = SS, Ref = Ref, Greedy = Greedy, P = P)]
    }

    # Maximization
    if(!is.null(solution)) {
      solution <- BiocGenerics::sort(solution)
      x2 <- x[data.table::inrange(x, IRanges::start(solution), IRanges::end(solution))]
      F2 <- ecdf(x2)
      S4Vectors::mcols(solution)$PSI <- mapply(function(i) mean(x2 %in% IRanges::start(solution[i]):IRanges::end(solution[i])), seq_along(solution))
      S4Vectors::mcols(solution)$beta <- mapply(function(i) {
        integral_result <- integrate(F2, lower = S4Vectors::mcols(solution[i])$SS - 1, upper = IRanges::end(solution[i]), subdivisions = 1000L)
        (integral_result$value - F2(S4Vectors::mcols(solution[i])$SS - 1) * (IRanges::end(solution[i]) - S4Vectors::mcols(solution[i])$SS + 1)) / ((IRanges::end(solution[i]) - S4Vectors::mcols(solution[i])$SS + 1) * (F2(IRanges::end(solution[i])) - F2(S4Vectors::mcols(solution[i])$SS - 1)))
      }, seq_along(solution))
    }
  } else {
    if(length(RefTSS) > 0 & is.element(L1[, pos], RefTSS)) {
      solution <- L1[, IRanges::IRanges(start = pos, end = pos, SS = pos, Ref = pos, Greedy = FALSE, P = 1, PSI = 1, beta = 1)]
    } else {
      solution <- L1[, IRanges::IRanges(start = pos, end = pos, SS = pos, Ref = NA, Greedy = FALSE, P = 1, PSI = 1, beta = 1)]
    }
  }

  # re-orientation
  if(!is.null(solution)) {
    solution <- data.table::as.data.table(solution)
    solution[, start := plyr::mapvalues(start, To, From, warn_missing = FALSE)]
    solution[, end := plyr::mapvalues(end, To, From, warn_missing = FALSE)]
    solution[, SS := plyr::mapvalues(SS, To, From, warn_missing = FALSE)]
    solution[, Ref := plyr::mapvalues(Ref, To, From, warn_missing = FALSE)]
    if(direction == "start") {
      solution <- solution[, IRanges::IRanges(start, end, SS = SS, Annotated = Ref, Greedy = Greedy, Percent = P, PSI = PSI, beta = beta)]
    } else {
      solution <- solution[, IRanges::IRanges(end, start, SS = SS, Annotated = Ref, Greedy = Greedy, Percent = P, PSI = PSI, beta = beta)]
    }
    solution <- BiocGenerics::sort(solution)
    if(length(RefTSS) != 0) RefTSS <- plyr::mapvalues(RefTSS, To, From, warn_missing = FALSE)
  }

  if(plots == 1) {
    par(mfrow = c(2, 1))
    L1[, plot(pos, CDR1, type = "s", ylim = c(0, 1))]
    L1[, lines(pos, CDR11, type = "s", col = 2)]
    L1[, lines(pos, CDR0, type = "s", col = "grey70")]
    L1[, plot(pos, d1, type = "l")]
    L1[, lines(pos, d0, type = "s", col = "grey70")]
    abline(h = ci, col = 3)
    abline(v = cd1[, Ref], lty = 2, col = "red")
    abline(v = cd1[, SS], lty = 2, col = "grey20")
  }

  if(plots == 2) {
    if(direction == "end") {
      y2 <- plyr::mapvalues(x, To, From, warn_missing = FALSE)
      y3 <- y2[y2 >= min(with(as.data.frame(solution), start)) & y2 <= max(with(as.data.frame(solution), end))]
      ggplot(data.table(x = min(with(as.data.frame(solution), start)):max(with(as.data.frame(solution), end)), y = ecdf(y3)(min(with(as.data.frame(solution), start)):max(with(as.data.frame(solution), end)))), aes(x, 1 - y)) +
        geom_step(colour = "grey") +
        scale_y_continuous(limits = c(0, 1)) +
        theme_classic() -> p1
    } else {
      y <- plyr::mapvalues(x, To, From, warn_missing = FALSE)
      ggplot(data.table(x = min(with(as.data.frame(solution), start)):max(with(as.data.frame(solution), end)), y = ecdf(y)(min(with(as.data.frame(solution), start)):max(with(as.data.frame(solution), end)))), aes(x, y)) +
        geom_step(colour = "grey") +
        scale_y_continuous(limits = c(0, 1)) +
        theme_classic() -> p1
    }
    if(!is.null(solution)) {
      p1 <- p1 + geom_vline(xintercept = with(as.data.frame(solution), SS), lty = 1)
      p1 <- p1 + geom_vline(xintercept = with(as.data.frame(solution), Annotated[!is.na(Annotated)]), lty = 2, colour = "red", lwd = 1.2)
    }
    print(p1)
  }
  return(solution)
}




#' @title InternalFunctions
#' @description Internal function for CliffSite

InflexionPoint <- function(x, y) {
  y2 <- S4Vectors::Rle(y)
  rV <- S4Vectors::runValue(y2)
  rL <- S4Vectors::runLength(y2)
  tp <- which(mapply(function(i) max(na.exclude(rV[(i-1):(i+1)])) == rV[i], seq_along(rV)))
  bp <- which(mapply(function(i) min(na.exclude(rV[(i-1):(i+1)])) == rV[i], seq_along(rV)))
  xy <- data.table::data.table(x, y, Index = rep(seq_along(rV), rL))
  tp <- do.call(rbind, lapply(tp, function(i) xy[Index == i, .(start = min(x), end = max(x), y = mean(y), pos = min(x), P = "T")]))
  bp <- do.call(rbind, lapply(bp, function(i) xy[Index == i, .(start = min(x), end = max(x), y = mean(y), pos = min(x), P = "B")]))
  rbind(bp, tp)[order(start)]
}

#' @title InternalFunctions
#' @description Internal function for CliffSite

DR <- function(x, minreads = 10, iteration = 100) {
  stopifnot(is(x, "data.table"))
  stopifnot(all(is.element(c("start", "end", "count", "index", "TSS"), colnames(x))))
  if(x[, max(count)] < minreads) return(NULL)
  x <- x[index %in% x[, .(count = max(count)), index][count >= minreads, index]]
  iso_g <- x[, unique(index)]
  mapply(function(i) {
    tot <- x[index == iso_g[i]]
    if(nrow(tot) != 0){
      tot <- tot[order(start, decreasing = T), ]
      INDEX <- matrix(1, nrow = nrow(tot), ncol = 1) # Igij
      DIST <- matrix(0, nrow = nrow(tot), ncol = 1)
      NM <- matrix(0, nrow = nrow(tot), ncol = 1)
      lexon <- tot[, end - start + 1] #lgj
      exon_c <- tot[, count]

      for(k in 1:nrow(tot)) {
        if(k == 1) {
          DIST[k, ] <- lexon[k] * INDEX[k, ] / 2
        }
        if(k > 1) {
          DIST[k, ] <- INDEX[k, ] * (lexon[1:k] %*% INDEX[1:k, ] - lexon[k] * INDEX[k, ] / 2)
        }
      }
      for(k in seq_len(ncol(DIST))) {
        #normalized distance
        DIST[, k] <- DIST[, k] / sum(lexon * INDEX[, k]) # sum(lexon*INDEX[, k] = Lgi; 得到dgij
      }

      ### EM
      alpha <- 0
      theta <- 1 / length(iso_g)
      for(iter in seq_len(iteration)) {
        ## E-Step
        if(is.na(alpha) | is.infinite(alpha)) break()
        EEE <- exp(- alpha * DIST)
        for(k in 1:nrow(tot)) {
          #print(k)
          eeek <- sum(INDEX[k, ] * theta * EEE[k, ]) # 公式2.5的分母部分第k个外显子在所有转录本表达*e的加和
          if(!is.na(eeek)) { # && ( eeek > 1e-12)
            NM[k, ] <- exon_c[k] * (INDEX[k, ] * theta * EEE[k, ]) / eeek #公式2.5 得到新的外显子count, 每个外显子在不同的转录本中count相等
          } else {
            NM[k, ] <- rep(0, length(iso_g))
          }
        }
        ## M-Step
        #if(is.finite(a) == 1){
        a <- alpha
        for(n in seq_len(iteration)) {
          # 迭代10次就可以让a值收敛（在supplementary里面公式（1）下面的一段话有说明）
          eee <- exp(- a * DIST)
          D0 <- lexon %*% (INDEX * eee) # 公式2.5下面的Gt的log函数的分子，在supplementary里面，即W0
          # %*% 是两个矩阵相乘
          D1 <- lexon %*% (INDEX * eee * DIST) # W1
          D2 <- lexon %*% (INDEX * eee * DIST^2) # W2
          TI <- colSums(INDEX * NM) # Tgi值，即每个转录本的read (位于supplementary公式（1）上的公式)
          First <- -sum(INDEX * NM * DIST) # 所有外显子的距离*read之和，（supplementary第二个公式的等式右边的第一项）
          Second <- 0
          for(ttt in 1:ncol(NM)) {
            # ttt是转录本的数量
            First <- First + TI[ttt] * (D1[ttt] / D0[ttt]) # supplementary第二个公式
            Second <- Second + TI[ttt] * ((D1[ttt]^2 - D0[ttt] * D2[ttt]) / D0[ttt]^2) # supplementary第三个公式
          }
          # print(paste0(n, "; ", First, "; ", Second))
          if(First == 0) break()
          a <- a - First / Second #公式（1）
        }
        if(identical(alpha, a)) break()
        alpha <- a
      }
      return(a)
    }
  }, seq_along(iso_g)) -> ATS_alpha
  names(ATS_alpha) <- iso_g

  if(length(unique(ATS_alpha)) != 0){ #第一次迭代失败，不开始第二次
    x <- x[index %in% names(ATS_alpha)]
    tot_g <- x[order(start, decreasing = T), ]
    iso_g <- tot_g[, unique(index)]
    INDEX <- matrix(0, nrow = nrow(tot_g), ncol = length(iso_g))
    DIST <- matrix(0, nrow = nrow(tot_g), ncol = length(iso_g))
    colnames(DIST) <- iso_g
    NM <- matrix(0, nrow = nrow(tot_g), ncol = length(iso_g))
    lexon <- tot_g[, end - start + 1] #lgj
    exon_c <- tot_g[, count] # Ngij
    for(k in 1:nrow(tot_g)) {
      #print(k)
      INDEX[k, ] <- as.numeric(iso_g <= tot_g[k, index]) # Igij
      if(k == 1) {
        DIST[k, ] <- lexon[k] * INDEX[k, ] / 2
      }
      if(k > 1) {
        DIST[k, ] <- INDEX[k, ] * (lexon[1:k] %*% INDEX[1:k, ] - lexon[k] * INDEX[k, ] / 2) #距离3‘越来越远
      }
    }
    for(k in 1:length(iso_g)) {
      #normalized distance
      DIST[, k] <- DIST[, k] / sum(lexon * INDEX[, k]) # sum(lexon * INDEX[, k] = Lgi; 得到dgij
    }

    alpha <- ATS_alpha[as.character(iso_g)]
    # theta <- rep(1 / length(iso_g), length(iso_g))
    theta <- prop.table(mapply(function(x) tot_g[index == x, max(count)], iso_g))

    for(iter in seq_len(iteration)) {
      ## E-Step
      if(anyNA(alpha) | any(is.infinite(alpha))) break() #一旦变成NA就不继续了，节约时间
      EEE <- exp(- alpha * DIST)
      for(k in 1:nrow(tot_g)) {
        # print(k)
        eeek <- sum(INDEX[k, ] * theta * EEE[k, ]) # 公式2.5的分母部分第k个外显子在所有转录本表达*e的加和
        if(!is.na(eeek) ) { # && ( eeek> 1e-12)
          NM[k, ] <- exon_c[k] * (INDEX[k, ] * theta * EEE[k, ]) / eeek
          # #公式2.5 得到新的外显子count, 每个外显子在不同的转录本中count相等
        } else {
          NM[k, ] <- rep(0, length(iso_g))
        }
      }
      ## M-Step
      a <- alpha
      for(n in seq_len(iteration)) { # 迭代10次就可以让a值收敛（在supplementary里面公式（1）下面的一段话有说明）
        # print(n)
        eee <- exp(- a * DIST)
        D0 <- lexon %*% (INDEX * eee) # 公式2.5下面的Gt的log函数的分子，在supplementary里面，即W0
        # %*% 是两个矩阵相乘
        D1 <- lexon %*% (INDEX * eee * DIST) # W1
        D2 <- lexon %*% (INDEX * eee * DIST^2) # W2
        TI <- colSums(INDEX * NM) # Tgi值，即每个转录本的read (位于supplementary公式（1）上的公式)
        First <- -sum(INDEX * NM * DIST) # 所有外显子的距离*read之和，（supplementary第二个公式的等式右边的第一项）
        Second <- 0
        for(ttt in 1:ncol(NM)) { # ttt是转录本的数量
          First <- First + TI[ttt] * (D1[ttt] / D0[ttt]) # supplementary第二个公式
          Second <- Second + TI[ttt] * ((D1[ttt]^2 - D0[ttt] * D2[ttt]) / D0[ttt]^2) # supplementary第三个公式
        }
        if(First == 0) break()
        a <- a - First / Second #公式（1）
      }
      if(identical(alpha, a)) break()
      alpha <- a
      ################## 求0gi （因为0gi的方程是关于a的，因此先求出a，再代入求0gi）
      EEE <- exp(- alpha * DIST) # 带入a的值得到，新的e对数项
      WWW <- lexon %*% (INDEX * EEE) # 进而得到wgh的值，在
      BBB <- (colSums(INDEX * NM) / sum(exon_c)) / WWW #2.6的分子
      theta <- BBB / sum(BBB) # 原文中公式2.6
    }
    if(!anyNA(alpha)) {
      # data.table(TSS = x[, levels(TSS)][iso_g], alpha1 = ATS_alpha[as.character(iso_g)], alpha2 = alpha, theta = as.numeric(theta))
      data.table(TSS = x[, levels(TSS)][iso_g], alpha = ATS_alpha[as.character(iso_g)], theta = as.numeric(theta))
    }
  }
}


#' @title InternalFunctions
#' @description Internal functions
myproportions <- function(x) {
  x[is.na(x)] <- 0
  x/sum(x)
}

#' @title InternalFunctions
#' @description Internal functions
loadbamregion <- function(bam, region, txdb = NULL,
                          exonIntersect = FALSE,
                          junctionIntersect = FALSE,
                          isSupplementaryAlignment = FALSE,
                          isSecondaryAlignment = FALSE,
                          isDuplicate = NA,
                          mapqFilter = 255,
                          ignore.strand = FALSE,
                          sep = "_") {
  options(warn = -1)
  stopifnot(all(file.exists(bam)))
  stopifnot(is.logical(exonIntersect))
  stopifnot(is.logical(junctionIntersect))
  if(exonIntersect) stopifnot(is(txdb, "TxDb"))
  if(junctionIntersect) stopifnot(is(txdb, "TxDb"))
  stopifnot(is.logical(isSupplementaryAlignment))
  stopifnot(is.logical(isSecondaryAlignment))
  stopifnot(is.logical(isDuplicate))
  stopifnot(is.logical(ignore.strand))
  stopifnot(is.numeric(mapqFilter))
  stopifnot(mapqFilter >= 0)
  stopifnot(is.character(sep))
  stopifnot(length(sep) == 1)

  param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = isSupplementaryAlignment,
                                                                 isPaired = TRUE,
                                                                 isFirstMateRead = TRUE,
                                                                 isDuplicate = isDuplicate,
                                                                 isNotPassingQualityControls = FALSE,
                                                                 isSecondaryAlignment = isSecondaryAlignment),
                                   tag = "CB", mapqFilter = mapqFilter)
  if(length(bam) == 1) {
    bami <- paste0(tempfile(), ".bam")
    sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, bam, gsub(":[-\\+]", "", as.character(region)))
    system(command = sh); system(command = paste("samtools index", bami))
    map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bami, param = param, use.names = TRUE, strandMode = 1))
    map0 <- map0[!is.na(S4Vectors::mcols(GenomicAlignments::first(map0))$CB)]
    if(length(map0) > 0) {
      S4Vectors::mcols(map0)$CB <- S4Vectors::mcols(GenomicAlignments::first(map0))$CB
    }
    system(command = paste0("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
  } else {
    map0 <- lapply(bam, function(x) {
      bami <- paste0(tempfile(), "_", basename(x))
      sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, x, gsub(":[-\\+]", "", as.character(region)))
      system(command = sh); system(command = paste("samtools index", bami))
      map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bami, param = param, use.names = TRUE, strandMode = 1))
      map0 <- map0[!is.na(S4Vectors::mcols(GenomicAlignments::first(map0))$CB)]
      if(length(map0) > 0) {
        S4Vectors::mcols(map0)$CB <- paste(S4Vectors::mcols(GenomicAlignments::first(map0))$CB, gsub(".bam", "", basename(x)), sep = sep)
      }
      system(command = paste("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
      return(map0)
    })
    map0 <- do.call(c, map0)
  }
  map0 <- map0[S4Vectors::queryHits(GenomicRanges::findOverlaps(GenomicRanges::GRanges(map0), region, ignore.strand = ignore.strand))]
  if(exonIntersect | junctionIntersect) {
    R1 <- GenomicAlignments::first(map0)
    GenomicRanges::strand(R1) <- GenomicRanges::strand(map0)
    R1_M <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = GenomicAlignments::cigar(R1), ops = "M", pos = GenomicAlignments::start(R1))
    names(R1_M) <- names(R1)
    R1_M <- unlist(R1_M)

    R2 <- GenomicAlignments::second(map0)
    GenomicRanges::strand(R2) <- GenomicRanges::strand(map0)
    R2_M <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = GenomicAlignments::cigar(R2), ops = "M", pos = GenomicAlignments::start(R2))
    names(R2_M) <- names(R2)
    R2_M <- unlist(R2_M)
  }

  if(exonIntersect) {
    R1Us <- names(R1_M)[S4Vectors::queryHits(IRanges::findOverlaps(R1_M, IRanges::ranges(GenomicFeatures::exonsByOverlaps(txdb, region)), type = "within"))]
    R2Us <- names(R2_M)[S4Vectors::queryHits(IRanges::findOverlaps(R2_M, IRanges::ranges(GenomicFeatures::exonsByOverlaps(txdb, region)), type = "within"))]
    map0 <- map0[names(map0) %in% union(R1Us, R2Us)]
  }

  if(junctionIntersect) {
    R1U1 <- names(R1_M)[S4Vectors::queryHits(IRanges::findOverlaps(R1_M, IRanges::ranges(GenomicFeatures::exonsByOverlaps(txdb, region)), type = "start"))]
    R1U2 <- names(R1_M)[S4Vectors::queryHits(IRanges::findOverlaps(R1_M, IRanges::ranges(GenomicFeatures::exonsByOverlaps(txdb, region)), type = "end"))]

    R2U1 <- names(R2_M)[S4Vectors::queryHits(IRanges::findOverlaps(R2_M, IRanges::ranges(GenomicFeatures::exonsByOverlaps(txdb, region)), type = "start"))]
    R2U2 <- names(R2_M)[S4Vectors::queryHits(IRanges::findOverlaps(R2_M, IRanges::ranges(GenomicFeatures::exonsByOverlaps(txdb, region)), type = "end"))]
    map0 <- map0[names(map0) %in% union(union(R1U1, R1U2), union(R2U1, R2U2))]
  }

  if(length(map0) == 0) {
    return(NULL)
  } else {
    return(map0)
  }
}

#' @title InternalFunctions
#' @description Internal functions
loadregionR1 <- function(bam, region,
                         isDuplicate = NA,
                         mapqFilter = 255,
                         ignore.strand = FALSE,
                         sep = "_") {
  options(warn = -1)
  stopifnot(all(file.exists(bam)))
  stopifnot(is.logical(ignore.strand))
  stopifnot(is.numeric(mapqFilter))
  stopifnot(mapqFilter >= 0)
  stopifnot(is.character(sep))
  stopifnot(length(sep) == 1)

  param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = FALSE,
                                                                 isPaired = TRUE,
                                                                 isFirstMateRead = TRUE,
                                                                 isDuplicate = isDuplicate,
                                                                 isNotPassingQualityControls = FALSE,
                                                                 isSecondaryAlignment = FALSE),
                                   tag = "CB", mapqFilter = mapqFilter)
  if(length(bam) == 1) {
    bami <- paste0(tempfile(), ".bam")
    sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, bam, gsub(":[-\\+]", "", as.character(region)))
    system(command = sh); system(command = paste("samtools index", bami))
    map0 <- suppressWarnings(GenomicAlignments::readGAlignments(bami, param = param, use.names = TRUE))
    map0 <- map0[!is.na(S4Vectors::mcols(map0)$CB)]
    if(length(map0) > 0) {
      S4Vectors::mcols(map0)$CB <- S4Vectors::mcols(map0)$CB
    }
    system(command = paste0("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
  } else {
    map0 <- lapply(bam, function(x) {
      bami <- paste0(tempfile(), "_", basename(x))
      sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, x, gsub(":[-\\+]", "", as.character(region)))
      system(command = sh); system(command = paste("samtools index", bami))
      map0 <- suppressWarnings(GenomicAlignments::readGAlignments(bami, param = param, use.names = TRUE))
      map0 <- map0[!is.na(S4Vectors::mcols(map0)$CB)]
      if(length(map0) > 0) {
        S4Vectors::mcols(map0)$CB <- paste(S4Vectors::mcols(map0)$CB, gsub(".bam", "", basename(x)), sep = sep)
      }
      system(command = paste("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
      return(map0)
    })
    map0 <- do.call(c, map0)
  }
  map0 <- map0[S4Vectors::queryHits(GenomicRanges::findOverlaps(GenomicRanges::GRanges(map0), region, ignore.strand = ignore.strand))]

  if(length(map0) == 0) {
    return(NULL)
  } else {
    return(map0)
  }
}

#' @title InternalFunctions
#' @description Internal functions
loadbamgene <- function(bam, gene, region = NULL, txdb = NULL, exonByGene = NULL,
                        exonIntersect = FALSE,
                        junctionIntersect = FALSE,
                        ignore.strand = FALSE,
                        isSupplementaryAlignment = FALSE,
                        isSecondaryAlignment = FALSE,
                        isDuplicate = NA,
                        mapqFilter = 255,
                        sep = "_") {
  options(warn = -1)
  stopifnot(length(gene) == 1)
  stopifnot(all(file.exists(bam)))
  stopifnot(is.logical(exonIntersect))
  stopifnot(is.logical(junctionIntersect))
  stopifnot(is.logical(isSupplementaryAlignment))
  stopifnot(is.logical(isSecondaryAlignment))
  stopifnot(is.logical(isDuplicate))
  stopifnot(is.logical(ignore.strand))
  stopifnot(is.numeric(mapqFilter))
  stopifnot(mapqFilter >= 0)
  stopifnot(is.character(sep))
  stopifnot(length(sep) == 1)
  if(is.null(region)) stopifnot(is(txdb, "TxDb"))
  if((exonIntersect | junctionIntersect) & is.null(exonByGene)) stopifnot(is(txdb, "TxDb"))

  if(is.null(region)) region <- GenomicFeatures::genes(txdb)[gene]
  param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = isSupplementaryAlignment,
                                                                 isPaired = TRUE,
                                                                 isDuplicate = isDuplicate,
                                                                 isNotPassingQualityControls = FALSE,
                                                                 isSecondaryAlignment = isSecondaryAlignment),
                                   tag = "CB", mapqFilter = mapqFilter)
  if(length(bam) == 1) {
    bami <- paste0(tempfile(), ".bam")
    sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, bam, gsub(":[-\\+]", "", as.character(region)))
    system(command = sh); system(command = paste("samtools index", bami))
    map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bami, param = param, use.names = TRUE, strandMode = 1))
    map0 <- map0[!is.na(S4Vectors::mcols(GenomicAlignments::first(map0))$CB)]
    if(length(map0) > 0) {
      S4Vectors::mcols(map0)$CB <- S4Vectors::mcols(GenomicAlignments::first(map0))$CB
    }
    system(command = paste0("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
  } else {
    map0 <- lapply(bam, function(x) {
      bami <- paste0(tempfile(), "_", basename(x))
      sh <- paste("samtools view -1h -q", mapqFilter, "-o", bami, x, gsub(":[-\\+]", "", as.character(region)))
      system(command = sh); system(command = paste("samtools index", bami))
      map0 <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bami, param = param, use.names = TRUE, strandMode = 1))
      map0 <- map0[!is.na(S4Vectors::mcols(GenomicAlignments::first(map0))$CB)]
      if(length(map0) > 0) {
        S4Vectors::mcols(map0)$CB <- paste(S4Vectors::mcols(GenomicAlignments::first(map0))$CB, gsub(".bam", "", basename(x)), sep = sep)
      }
      system(command = paste("rm ", bami)); system(command = paste0("rm ", bami, ".bai"))
      return(map0)
    })
    map0 <- do.call(c, map0)
  }
  map0 <- map0[S4Vectors::queryHits(GenomicRanges::findOverlaps(GenomicRanges::GRanges(map0), region, ignore.strand = ignore.strand))]
  if(exonIntersect | junctionIntersect) {
    if(!is.null(exonByGene)) {
      if(is(exonByGene, "CompressedGRangesList")) {
        if(is.element(gene, names(exonByGene))) {
          exonByGene <- exonByGene[[gene]]
        } else {
          exonByGene <- GenomicFeatures::exonsBy(txdb, "gene")[[gene]]
        }
      } else {
        exonByGene <- GenomicFeatures::exonsBy(txdb, "gene")[[gene]]
      }
    } else {
      exonByGene <- GenomicFeatures::exonsBy(txdb, "gene")[[gene]]
    }

    R1 <- GenomicAlignments::first(map0)
    GenomicRanges::strand(R1) <- GenomicRanges::strand(map0)
    R1_M <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = GenomicAlignments::cigar(R1), ops = "M", pos = GenomicAlignments::start(R1))
    names(R1_M) <- names(R1)
    R1_M <- unlist(R1_M)

    R2 <- GenomicAlignments::second(map0)
    GenomicRanges::strand(R2) <- GenomicRanges::strand(map0)
    R2_M <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = GenomicAlignments::cigar(R2), ops = "M", pos = GenomicAlignments::start(R2))
    names(R2_M) <- names(R2)
    R2_M <- unlist(R2_M)
  }
  if(exonIntersect) {
    R1Us <- names(R1_M)[S4Vectors::queryHits(IRanges::findOverlaps(R1_M, IRanges::ranges(exonByGene), type = "within"))]
    R2Us <- names(R2_M)[S4Vectors::queryHits(IRanges::findOverlaps(R2_M, IRanges::ranges(exonByGene), type = "within"))]
    map0 <- map0[names(map0) %in% union(R1Us, R2Us)]
  }
  if(junctionIntersect) {
    R1U1 <- names(R1_M)[S4Vectors::queryHits(IRanges::findOverlaps(R1_M, IRanges::ranges(exonByGene), type = "start"))]
    R1U2 <- names(R1_M)[S4Vectors::queryHits(IRanges::findOverlaps(R1_M, IRanges::ranges(exonByGene), type = "end"))]

    R2U1 <- names(R2_M)[S4Vectors::queryHits(IRanges::findOverlaps(R2_M, IRanges::ranges(exonByGene), type = "start"))]
    R2U2 <- names(R2_M)[S4Vectors::queryHits(IRanges::findOverlaps(R2_M, IRanges::ranges(exonByGene), type = "end"))]
    map0 <- map0[names(map0) %in% union(union(R1U1, R1U2), union(R2U1, R2U2))]
  }
  if(length(map0) == 0) {
    return(NULL)
  } else {
    return(map0)
  }
}

#' @title InternalFunctions
#' @description Internal functions

CliffSite2 <- function(x, SS, direction = "start", minreads = 100, max.degradation.length = NULL, min.Nodegradation.length = NULL) {
  stopifnot(is.element(direction, c("start", "end")))
  if(length(x) < minreads) return(NULL)
  From <- sort(unique(c(x, start(SS), end(SS), S4Vectors::mcols(SS)$SS)))
  # To <- seq_along(From)
  To <- From
  x <- plyr::mapvalues(x, From, To, warn_missing = FALSE)
  solution <- IRanges::IRanges(plyr::mapvalues(start(SS), From, To, warn_missing = FALSE),
                               plyr::mapvalues(end(SS), From, To, warn_missing = FALSE),
                               SS = plyr::mapvalues(S4Vectors::mcols(SS)$SS, From, To, warn_missing = FALSE))
  if(direction == "end") {
    z <- c(min(x) - 10, max(x) + 10)
    y <- x
    x <- min(z) + max(z) - x
    solution <- IRanges::IRanges(min(z) + max(z) - end(solution),
                                 min(z) + max(z) - start(solution),
                                 SS = min(z) + max(z) - S4Vectors::mcols(solution)$SS)
  }
  F1 <- ecdf(x)
  # Maximization
  if(!is.null(solution)) {
    solution <- BiocGenerics::sort(solution)
    mcols(solution)$PSI <- mapply(function(i) F1(with(as.data.frame(solution[i]), end)) - F1(with(as.data.frame(solution[i]), start)), seq_along(solution))
    if(!is.null(max.degradation.length) | !is.null(min.Nodegradation.length)) {
      x2 <- lapply(function(i) {
        xi <- x[data.table::between(x, with(as.data.frame(solution[i]), SS), with(as.data.frame(solution[i]), end))]
        if(!is.null(max.degradation.length)) {
          stopifnot(is.integer(as.integer(max.degradation.length)))
          xi <- xi[xi <= min(xi) + max.degradation.length]
        }
        if(!is.null(min.Nodegradation.length)) {
          si <- IRanges::gaps(IRanges::reduce(IRanges::IRanges(xi)))
          if(max(width(si)) > min.Nodegradation.length) {
            xi <- xi[xi <= min(start(si[width(si) > min.Nodegradation.length]))]
          }
        }
        xi
      }, X = seq_along(solution))
      x2 <- do.call(c, x2)
      F2 <- ecdf(x2)
      DRT <- lapply(function(i) {
        xi <- x[data.table::between(x, with(as.data.frame(solution[i]), SS), with(as.data.frame(solution[i]), end))]
        if(length(xi) == 0) return(NULL)
        if(length(unique(xi)) == 1) return(NULL)
        data.table(start = sort(unique(xi)), end = sort(unique(xi)), count = F2(sort(unique(xi))) * length(x2), index = with(as.data.frame(solution[i]), paste0(start, "-", end)))
      }, X = seq_along(solution))
      DRT <- do.call(rbind, DRT)
      DRR <- if(is.null(DRT)) NULL else DR(x = DRT[, .SD[count != max(count)], index])
    } else {
      x2 <- x
      DRT <- lapply(function(i) {
        xi <- x[data.table::between(x, with(as.data.frame(solution[i]), SS), with(as.data.frame(solution[i]), end))]
        if(length(xi) == 0) return(NULL)
        if(length(unique(xi)) == 1) return(NULL)
        data.table(start = sort(unique(xi)), end = sort(unique(xi)), count = F1(sort(unique(xi))) * length(x2), index = with(as.data.frame(solution[i]), paste0(start, "-", end)))
      }, X = seq_along(solution))
      DRT <- do.call(rbind, DRT)
      DRR <- if(is.null(DRT)) NULL else DR(x = DRT)
    }

    mi <- lapply(function(i) {
      xi <- x2[data.table::between(x2, with(as.data.frame(solution[i]), SS), with(as.data.frame(solution[i]), end))]
      if(length(xi) == 0) {
        return(data.table(AUC = NA, MaxDist = NA, AboveRandom = NA, ApicesX = NA, ApicesY = NA))
      }
      if(length(unique(xi)) == 1) {
        return(data.table(AUC = 1, MaxDist = 1, AboveRandom = 1, ApicesX = 0, ApicesY = 1))
      }
      Fi <- ecdf(xi)
      T0 <- data.table(x = min(xi):max(xi), E = seq(0, 1, length.out = length(min(xi):max(xi))), O = Fi(min(xi):max(xi)))
      data.table(AUC = T0[, sum(O) / .N],
                 MaxDist = T0[, max(O - E)],
                 AboveRandom = T0[, mean(O >= E)],
                 ApicesX = T0[, which.max(O - E) / .N],
                 ApicesY = T0[which.max(O - E), O])
    }, X = seq_along(solution))
    mi <- do.call(rbind, mi)
    mcols(solution) <- cbind(mcols(solution), mi)

    if(!is.null(DRR)) {
      setkey(DRR, ATS)
      mcols(solution)$alpha1 <- DRR[with(as.data.frame(solution), paste0(start, "-", end)), alpha1]
      mcols(solution)$alpha2 <- DRR[with(as.data.frame(solution), paste0(start, "-", end)), alpha2]
      mcols(solution)$theta <- DRR[with(as.data.frame(solution), paste0(start, "-", end)), theta]
    } else {
      mcols(solution)$alpha1 <- NA
      mcols(solution)$alpha2 <- NA
      mcols(solution)$theta <- NA
    }
  }

  # re-orientation
  if(!is.null(solution)) {
    if(direction == "end") {
      solution <- as.data.table(solution)
      solution[, start := plyr::mapvalues(start, min(z):max(z), max(z):min(z), warn_missing = FALSE)]
      solution[, end := plyr::mapvalues(end, min(z):max(z), max(z):min(z), warn_missing = FALSE)]
      solution[, SS := plyr::mapvalues(SS, min(z):max(z), max(z):min(z), warn_missing = FALSE)]
      solution <- solution[, IRanges(start = end, end = start, SS = SS, PSI = PSI, AUC = AUC, MaxDist = MaxDist, AboveRandom = AboveRandom, ApicesX = ApicesX, ApicesY = ApicesY, alpha1 = alpha1, alpha2 = alpha2, theta = theta)]
      solution <- BiocGenerics::sort(solution)
    }
  }

  #  Window-switch
  if(length(solution) > 0) {
    solution <- as.data.table(solution)
    solution[, start := plyr::mapvalues(start, To, From, warn_missing = FALSE)]
    solution[, end := plyr::mapvalues(end, To, From, warn_missing = FALSE)]
    solution[, SS := plyr::mapvalues(SS, To, From, warn_missing = FALSE)]
    solution <- solution[, IRanges(start, end, SS = SS, PSI = PSI, AUC = AUC, MaxDist = MaxDist, AboveRandom = AboveRandom, ApicesX = ApicesX, ApicesY = ApicesY, alpha1 = alpha1, alpha2 = alpha2, theta = theta)]
    solution <- BiocGenerics::sort(solution)
  }
  SS <- sort(SS)
  stopifnot(identical(as.character(solution), as.character(SS)))
  return(solution)
}

#' @title InternalFunctions
#' @description Internal functions
determineWindowLength <- function (x, w) {
  if (!is.numeric(x) | !is.numeric(w)) {
    stop("Arguments 'x' and 'w' must be numeric")
  }
  l = length(x)
  w.orig = w
  if (w > 0 && w < 1) {
    w = w * l
  }
  ret = as.integer(max(abs(w[1]), 1))
  if (length(x) <= ret | ret < 1) {
    stop(paste("Resultant window length is out of range (", ret, "), must be >= 1", sep = ""), call. = FALSE)
  }
  return(ret)
}

#' @title InternalFunctions
#' @description Internal functions
normalize <- function (x) {
  if (!is.numeric(x)) {
    stop("argument 'x' must be numeric")
  }
  if (length(x) == 1) {
    return(1)
  }
  total = sum(x)
  if (total == 0) {
    stop("argument 'x' must not sum to zero")
  }
  return(x/total)
}

#' @title InternalFunctions
#' @description Internal functions
smthgaussian <- function (x, window = 0.1, alpha = 2.5, tails = TRUE) {
  if (!is.numeric(x) | !is.numeric(alpha)) {
    stop("argument 'x' and 'alpha' must be numeric", call. = FALSE)
  }
  windowLength = determineWindowLength(x, window)
  makeWindow = function(w, a) {
    hw = abs(w/2)
    e = exp(1)
    a = abs(a)
    ret = sapply(c(0:(w - 1)), function(x) {
      n = x - as.integer(hw)
      k = -0.5 * (a * n/hw)^2
      e^k
    })
    ret
  }
  w = makeWindow(windowLength, alpha[1])
  sizeW = length(w)
  sizeD = length(x)
  w = normalize(w)
  hkwL = as.integer(sizeW/2)
  hkwR = sizeW - hkwL
  ret = sapply(c(1:sizeD), function(i) {
    ix.d = c((i - hkwL):(i + hkwR - 1))
    ix.w = which(ix.d %in% 1:sizeD)
    ix.d = ix.d[ix.w]
    W.nm = if (length(ix.w) != sizeW) {
      normalize(w[ix.w])
    }
    else {
      w
    }
    D.nm = x[ix.d]
    as.numeric(D.nm %*% W.nm)
  })
  if (!tails) {
    ret[c(1:hkwL, (sizeD - hkwR + 1):sizeD)] = NA
  }
  ret
}
