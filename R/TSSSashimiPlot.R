#' @title Sashimi
#' @name Sashimi
#' @description A function for TSS Sashimi plot.
#'
#' @importFrom rtracklayer readGFFAsGRanges
#' @importFrom S4Vectors queryHits mcols
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom GenomicAlignments readGAlignmentPairs second first findOverlaps coverage
#' @importFrom IRanges start ranges end
#' @importFrom data.table data.table
#' @importFrom GenomicRanges seqnames start end reduce gaps GRangesList findOverlaps
#' @importFrom grDevices colorRampPalette col2rgb
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom IRanges IRanges start end findOverlaps ranges
#' @importFrom plyr mapvalues
#' @importFrom ggbio autoplot tracks
#' @importFrom patchwork plot_layout
#' @import ggplot2
#' @import patchwork
#'
#' @param object An scATSDataSet object.
#' @param reads The mapped reads, imported by loadbam function.
#' @param gtf The file path pf GTF file, or GRanges of GTF file imported by readGFFAsGRanges of rtracklayer package.
#' @param gene The target gene.
#' @param which A GRanges, or any object that can be coerced to a GRanges, or missing object.
#' @param Cells Which cells will by remained.
#' @param groupBy,groupLevels Group the cells, and which groups will be used for Sashimi.
#' @param transcripts Which transcripts will by remained in the reference panel.
#' @param TSS The loci of TSSs to be indicated in in Sashimi.
#' @param fill.color The color(s) of Sashimi of each group.
#' @param line.type,line.color The line type and line colour of indicated TSS.
#' @param base_size base font size of ggplot theme, given in pts.
#' @param rel_height relative height of reference panel.
#' @param free_y Free the y dimension of grouped Sashimi?
#' @param adjust,bw the smoothing bandwidth to be used of density ridges.
#'
#' @export
#' @concept visualization

Sashimi <- function(object, bam, sep = "_", gtf, txdb = NULL, gene = NULL, Cells = NULL, transcripts = NULL, xlimit = NULL, MinDepth = 20,
                    CellBarcodeName = NULL, groupBy = NULL, groups = NULL, groupLevels = NULL, TSS = NULL,
                    exonIntersect = TRUE, junctionIntersect = FALSE, ignore.strand = FALSE, mapqFilter = 255,
                    isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = NA,
                    fill.color = NULL, line.type = NULL, line.color = NULL, ColourTxName = "Red", base_size = 15, rel_height = NULL, free_y = TRUE,
                    Group.name.angle = NULL) {
  options(warn = -1)
  stopifnot(is(object, "scATSDataSet"))
  stopifnot(all(file.exists(bam)))
  stopifnot(all(file.exists(paste0(bam, ".bai"))))
  stopifnot(!is.null(gene))
  if(length(bam) > 1) {
    stopifnot(!is.null(sep))
    stopifnot(!is.na(sep))
    stopifnot(is.character(sep))
    stopifnot(length(sep) == 1)
    SampleName <- sort(unique(mapply(function(x) x[2], strsplit(colnames(object), split = sep))))
    stopifnot(identical(SampleName, sort(gsub(".bam", "", basename(bam)))))
  }
  if(is.character(gtf)) {
    stopifnot(file.exists(gtf))
    if(is.null(txdb) | !is(txdb, "TxDb")) {
      txdb <- suppressMessages(suppressWarnings(GenomicFeatures::makeTxDbFromGFF(gtf)))
    }
    gtf <- rtracklayer::readGFFAsGRanges(gtf)
    stopifnot(all(c("gene_name", "type") %in% names(S4Vectors::mcols(gtf))))
  } else {
    stopifnot(is(gtf, "GRanges"))
    stopifnot(all(c("gene_name", "type") %in% names(S4Vectors::mcols(gtf))))
    stopifnot(is(txdb, "TxDb"))
  }
  if(is.null(gene)) {
    stop("There is no gene for ATS identification!")
  } else {
    stopifnot(is.element(gene, SummarizedExperiment::rowData(object)$gene_name))
    stopifnot(is.element(gene, gtf$gene_name))
  }
  if(!is.null(Cells)) {
    Cells <- Cells[Cells %in% colnames(object)]
    stopifnot(length(Cells) != 0)
    object <- object[, Cells]
  } else {
    Cells <- colnames(object)
  }

  if(!is.null(xlimit)) {
    stopifnot(length(xlimit) == 2)
    stopifnot(xlimit[2] > xlimit[1])
  }

  stopifnot(MinDepth >= 0)

  if(!is.null(groupBy)) {
    stopifnot(is.element(groupBy, colnames(SummarizedExperiment::colData(object))))
    if(!is.null(CellBarcodeName)) {
      stopifnot(is.element(CellBarcodeName, colnames(SummarizedExperiment::colData(object))))
    } else {
      CellBarcodeName <- "arn"
    }
  } else {
    if(!is.null(CellBarcodeName)) {
      stopifnot(is.element(CellBarcodeName, colnames(SummarizedExperiment::colData(object))))
    } else {
      CellBarcodeName <- "arn"
    }
  }
  if(!is.null(groups)) stopifnot(all(is.element(groups, unique(subset(SummarizedExperiment::colData(object), select = groupBy)[[1]]))))
  if(!is.null(groupBy) & is.null(groups)) groups <- sort(unique(subset(SummarizedExperiment::colData(object), select = groupBy)[[1]]))
  if(!is.null(groupLevels)) stopifnot(all(is.element(groups, groups)))
  if(is.null(groupLevels)) groupLevels <- groups
  groups <- groups[groups %in% groupLevels]

  stopifnot(is.logical(exonIntersect))
  stopifnot(is.logical(junctionIntersect))
  stopifnot(is.logical(isSupplementaryAlignment))
  stopifnot(is.logical(isSecondaryAlignment))
  stopifnot(is.logical(isDuplicate))
  stopifnot(is.logical(ignore.strand))
  stopifnot(is.numeric(mapqFilter))
  stopifnot(mapqFilter >= 0)

  if(!is.null(line.type)) {
    stopifnot(is.numeric(line.type))
    stopifnot(length(line.type) == 1)
  } else {
    line.type <- 4
  }

  if(!is.null(line.color)) {
    stopifnot(is.character(line.color))
    stopifnot(length(line.color) == 1)
  } else {
    line.color <- "black"
  }

  if(!is.null(fill.color)) {
    stopifnot(is.character(fill.color))
  }

  if(!is.null(rel_height)) {
    stopifnot(length(rel_height) == 1)
    stopifnot(is.numeric(rel_height))
    stopifnot(rel_height > 0)
  }

  stopifnot(is.logical(free_y))

  if(!is.null(Group.name.angle)) {
    stopifnot(length(Group.name.angle) == 1)
    stopifnot(is.numeric(Group.name.angle))
  }

  EBG <- GenomicFeatures::exonsBy(txdb, by = "gene")

  paras <- subset(object@rowRanges, gene_name == gene)
  cellmeta <- data.table::as.data.table(SummarizedExperiment::colData(object), keep.rownames = "arn")[, c(groupBy, CellBarcodeName), with = F]
  eval(parse(text = paste0("data.table::setkey(cellmeta, ", CellBarcodeName, ")")))

  if(is.logical(TSS)) {
    if(TSS) {
      TSS <- GenomicRanges::start(paras)
    } else {
      TSS <- NULL
    }
  }

  region <- gtf[with(data.frame(S4Vectors::mcols(gtf)), type == "gene" & gene_name == gene)]
  reads <- loadbamgene(bam = bam,
                       gene = S4Vectors::mcols(region)$gene_id,
                       region = region,
                       txdb = txdb,
                       exonByGene = EBG,
                       exonIntersect = exonIntersect,
                       junctionIntersect = junctionIntersect,
                       ignore.strand = FALSE,
                       isSupplementaryAlignment = isSupplementaryAlignment,
                       isSecondaryAlignment = isSecondaryAlignment,
                       isDuplicate = isDuplicate,
                       mapqFilter = mapqFilter,
                       sep = sep)
  reads <- reads[S4Vectors::mcols(reads)$CB %in% Cells]

  if(!is.null(groupBy)) {
    cov <- split(colnames(object), SummarizedExperiment::colData(object)[, groupBy])
    reads <- lapply(cov, function(x) {
      reads[S4Vectors::mcols(reads)$CB %in% as.character(x)]
    })
  } else {
    reads <- list("Depth" = reads)
  }
  reads <- lapply(reads, function(x) {
    if(S4Vectors::runValue(GenomicAlignments::strand(region)) == "+") {
      Pos <- data.table::data.table(pos = IRanges::start(IRanges::ranges(x)), CB = S4Vectors::mcols(x)$CB)
      cov <- GenomicAlignments::coverage(GenomicAlignments::first(x))[[unique(as.character(GenomeInfoDb::seqnames(region)))]]
    } else {
      Pos <- data.table::data.table(pos = IRanges::end(IRanges::ranges(x)), CB = S4Vectors::mcols(x)$CB)
      cov <- GenomicAlignments::coverage(GenomicAlignments::first(x))[[unique(as.character(GenomeInfoDb::seqnames(region)))]]
    }
    if(nrow(Pos) > 0) {
      cov <- data.table::data.table(pos = Pos[, min(pos)]:Pos[, max(pos)],
                                    D = as.numeric(cov[Pos[, min(pos)]:Pos[, max(pos)]]))
    } else {
      Pos <- data.table::data.table(pos = NA, CB = NA)
      cov <- data.table::data.table(pos = NA, D = NA)
    }
    return(list(cov, Pos))
  })
  Pos <- lapply(reads, function(x) x[[2]])
  cov <- lapply(reads, function(x) x[[1]])
  Pos <- data.table::data.table(groupBy = rep(names(Pos), mapply(nrow, Pos)), do.call(rbind, Pos))
  cov <- data.table::data.table(groupBy = rep(names(cov), mapply(nrow, cov)), do.call(rbind, cov))
  if(!is.null(groupLevels)) cov[, groupBy := factor(groupBy, levels = groupLevels)]
  if(!is.null(groupLevels)) Pos[, groupBy := factor(groupBy, levels = groupLevels)]
  cov <- cov[groupBy %in% cov[, max(D), groupBy][V1 >= MinDepth, groupBy]]

  if(is.null(fill.color)) {
    fill.color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))(length(cov[, unique(groupBy)]))
  } else {
    if(length(fill.color) == 1) {
      if(fill.color %in% row.names(RColorBrewer::brewer.pal.info)) {
        fill.color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info[fill.color, "maxcolors"], fill.color))(length(cov[, unique(groupBy)]))
      } else {
        e <- tryCatch(col2rgb(fill.color), error = function(e) e)
        if(is(e, "simpleError")) {
          message("invalid color name, using grey to replace.")
          fill.color <- rep("grey", length(cov[, unique(groupBy)]))
        } else {
          fill.color <- rep(fill.color, length(cov[, unique(groupBy)]))
        }
      }
    } else {
      if(length(fill.color) != length(cov[, unique(groupBy)])) {
        message("invalid length of fill.color, using RColorBrewer Set1 to replace.")
        fill.color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))(length(cov[, unique(groupBy)]))
      }
    }
  }

  if(is.null(xlimit)) {
    xlimit <- c(min(c(GenomicRanges::start(region), TSS)), max(c(GenomicRanges::end(region), TSS)))
  }

  gi <- biovizBase::crunch(txdb, region)
  gi <- gi[(GenomicRanges::start(gi) >= min(xlimit) & GenomicRanges::start(gi) <= max(xlimit)) | (GenomicRanges::end(gi) >= min(xlimit) & GenomicRanges::end(gi) <= max(xlimit))]
  GenomicRanges::start(gi[GenomicRanges::start(gi) < min(xlimit)]) <- min(xlimit)
  GenomicRanges::end(gi[GenomicRanges::end(gi) > max(xlimit)]) <- max(xlimit)
  gi <- GenomicRanges::GRangesList(split(gi, gi$tx_name))
  names(gi) <- plyr::mapvalues(names(gi), with(as.data.frame(gtf), transcript_id[type == "transcript"]), with(as.data.frame(gtf), transcript_name[type == "transcript"]), warn_missing = F)
  if(!is.null(transcripts)) {
    if(any(transcripts %in% names(gi))) {
      gi <- gi[names(gi) %in% transcripts]
    }
  }
  gi <- gi[mapply(length, gi) > 0]
  id <- mapply(function(x) paste(sort(c(start(x), end(x))), collapse = ","), gi)
  if(any(duplicated(id))) {
    id <- split(names(gi), id); names(id) <- NULL
    gi <- lapply(id, function(x) {
      gi[[x[1]]]
    })
    names(gi) <- mapply(function(x) paste(x, collapse = "; "), id)
  }
  gi <- GenomicRanges::GRangesList(gi)
  p0 <- ggbio::autoplot(gi)@ggplot + theme_void()

  cov <- cov[pos %between% ggplot2::layer_scales(p0)$x$range$range]
  if(is.null(groupBy)) {
    if(is.null(Group.name.angle)) Group.name.angle <- 90
    ggplot(cov) +
      geom_ribbon(aes(x = pos, ymax = D, ymin = 0, fill = groupBy)) +
      scale_x_continuous(limits = ggplot2::layer_scales(p0)$x$range$range) +
      scale_y_continuous(n.breaks = 3) +
      scale_fill_manual(values = fill.color) +
      theme_light(base_size = base_size) +
      theme(legend.position = "none",
            axis.title = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.line.y = element_line(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.y.left = element_text(angle = Group.name.angle, colour = "black"),
            strip.placement = "outside",
            strip.background = element_rect(fill = NA)) -> p1
  } else {
    if(is.null(Group.name.angle)) Group.name.angle <- 0
    if(free_y) {
      ggplot(cov) +
        geom_ribbon(aes(x = pos, ymax = D, ymin = 0, fill = groupBy)) +
        facet_wrap(~ groupBy, ncol = 1, strip.position = "left", scales = "free_y") +
        scale_x_continuous(limits = ggplot2::layer_scales(p0)$x$range$range) +
        scale_y_continuous(n.breaks = 3) +
        scale_fill_manual(values = fill.color) +
        theme_light(base_size = base_size) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(),
              axis.line.y = element_line(),
              axis.ticks.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.y.left = element_text(angle = Group.name.angle, colour = "black"),
              strip.placement = "outside",
              strip.background = element_rect(fill = NA)) -> p1
    } else {
      ggplot(cov) +
        geom_ribbon(aes(x = pos, ymax = D, ymin = 0, fill = groupBy)) +
        facet_wrap(~ groupBy, ncol = 1, strip.position = "left") +
        scale_x_continuous(limits = ggplot2::layer_scales(p0)$x$range$range) +
        scale_y_continuous(n.breaks = 3) +
        scale_fill_manual(values = fill.color) +
        theme_light(base_size = base_size) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(),
              axis.line.y = element_line(),
              axis.ticks.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.y.left = element_text(angle = Group.name.angle, colour = "black"),
              strip.placement = "outside",
              strip.background = element_rect(fill = NA)) -> p1
    }
  }

  if(!is.null(TSS)) {
    p0 <- p0 + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
    p1 <- p1 + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
  }

  if(is.null(rel_height)) {
    rh <- c(cov[, length(unique(groupBy))], max(0.5, length(gi)/10))
    rh <- rh / max(rh)
    rh[rh < 0.2] <- 0.2
    # ggbio::tracks(p1, p4@ggplot, heights = rh, padding = unit(0, "lines")) + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
    p1 / p0 + patchwork::plot_layout(ncol = 1, heights = rh)
  } else {
    # ggbio::tracks(p1, p4@ggplot, heights = c(1, rel_height), padding = unit(0, "lines")) + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
    p1 / p0 + patchwork::plot_layout(ncol = 1, heights = c(1, rel_height))
  }
}


#' @title Sashimi
#' @name Sashimi
#'
#' @export
#' @concept visualization

StartSite <- function(object, bam, sep = "_", gtf, txdb = NULL, gene = NULL, Cells = NULL, transcripts = NULL, xlimit = NULL, MinDepth = 20,
                      CellBarcodeName = NULL, groupBy = NULL, groups = NULL, groupLevels = NULL, TSS = NULL,
                      exonIntersect = TRUE, junctionIntersect = FALSE, ignore.strand = FALSE, mapqFilter = 255,
                      isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = NA,
                      point.color = "#FEB24C", point.size = 1, point.shape = 16, line.type = NULL, line.color = NULL, ColourTxName = "Red",
                      base_size = 15, rel_height = NULL, free_y = FALSE, Group.name.angle = NULL) {
  options(warn = -1)
  stopifnot(is(object, "scATSDataSet"))
  stopifnot(all(file.exists(bam)))
  stopifnot(all(file.exists(paste0(bam, ".bai"))))
  stopifnot(!is.null(gene))
  if(length(bam) > 1) {
    stopifnot(!is.null(sep))
    stopifnot(!is.na(sep))
    stopifnot(is.character(sep))
    stopifnot(length(sep) == 1)
    SampleName <- sort(unique(mapply(function(x) x[2], strsplit(colnames(object), split = sep))))
    stopifnot(identical(SampleName, sort(gsub(".bam", "", basename(bam)))))
  }
  if(is.character(gtf)) {
    stopifnot(file.exists(gtf))
    if(is.null(txdb) | !is(txdb, "TxDb")) {
      txdb <- suppressMessages(suppressWarnings(GenomicFeatures::makeTxDbFromGFF(gtf)))
    }
    gtf <- rtracklayer::readGFFAsGRanges(gtf)
    stopifnot(all(c("gene_name", "type") %in% names(S4Vectors::mcols(gtf))))
  } else {
    stopifnot(is(gtf, "GRanges"))
    stopifnot(all(c("gene_name", "type") %in% names(S4Vectors::mcols(gtf))))
    stopifnot(is(txdb, "TxDb"))
  }
  if(is.null(gene)) {
    stop("There is no gene for ATS identification!")
  } else {
    stopifnot(is.element(gene, SummarizedExperiment::rowData(object)$gene_name))
    stopifnot(is.element(gene, gtf$gene_name))
  }
  if(!is.null(Cells)) {
    Cells <- Cells[Cells %in% colnames(object)]
    stopifnot(length(Cells) != 0)
    object <- object[, Cells]
  } else {
    Cells <- colnames(object)
  }

  if(!is.null(xlimit)) {
    stopifnot(length(xlimit) == 2)
    stopifnot(xlimit[2] > xlimit[1])
  }

  stopifnot(MinDepth >= 0)

  if(!is.null(groupBy)) {
    stopifnot(is.element(groupBy, colnames(SummarizedExperiment::colData(object))))
    if(!is.null(CellBarcodeName)) {
      stopifnot(is.element(CellBarcodeName, colnames(SummarizedExperiment::colData(object))))
    } else {
      CellBarcodeName <- "arn"
    }
  } else {
    if(!is.null(CellBarcodeName)) {
      stopifnot(is.element(CellBarcodeName, colnames(SummarizedExperiment::colData(object))))
    } else {
      CellBarcodeName <- "arn"
    }
  }
  if(!is.null(groups)) stopifnot(all(is.element(groups, unique(subset(SummarizedExperiment::colData(object), select = groupBy)[[1]]))))
  if(!is.null(groupBy) & is.null(groups)) groups <- sort(unique(subset(SummarizedExperiment::colData(object), select = groupBy)[[1]]))
  if(!is.null(groupLevels)) stopifnot(all(is.element(groups, groups)))
  if(is.null(groupLevels)) groupLevels <- groups
  groups <- groups[groups %in% groupLevels]

  stopifnot(is.logical(exonIntersect))
  stopifnot(is.logical(junctionIntersect))
  stopifnot(is.logical(isSupplementaryAlignment))
  stopifnot(is.logical(isSecondaryAlignment))
  stopifnot(is.logical(isDuplicate))
  stopifnot(is.logical(ignore.strand))
  stopifnot(is.numeric(mapqFilter))
  stopifnot(mapqFilter >= 0)

  if(!is.null(line.type)) {
    stopifnot(is.numeric(line.type))
    stopifnot(length(line.type) == 1)
  } else {
    line.type <- 4
  }

  if(!is.null(line.color)) {
    stopifnot(is.character(line.color))
    stopifnot(length(line.color) == 1)
  } else {
    line.color <- "black"
  }

  if(!is.null(point.color)) {
    stopifnot(is.character(point.color))
  }
  stopifnot(point.size > 0)
  stopifnot(point.shape >= 0 & point.shape <= 25)

  if(!is.null(rel_height)) {
    stopifnot(length(rel_height) == 1)
    stopifnot(is.numeric(rel_height))
    stopifnot(rel_height > 0)
  }

  stopifnot(is.logical(free_y))

  if(!is.null(Group.name.angle)) {
    stopifnot(length(Group.name.angle) == 1)
    stopifnot(is.numeric(Group.name.angle))
  }

  EBG <- GenomicFeatures::exonsBy(txdb, by = "gene")

  paras <- subset(object@rowRanges, gene_name == gene)
  cellmeta <- data.table::as.data.table(SummarizedExperiment::colData(object), keep.rownames = "arn")[, c(groupBy, CellBarcodeName), with = F]
  eval(parse(text = paste0("data.table::setkey(cellmeta, ", CellBarcodeName, ")")))

  if(is.logical(TSS)) {
    if(TSS) {
      TSS <- GenomicRanges::start(paras)
    } else {
      TSS <- NULL
    }
  }

  region <- gtf[with(data.frame(S4Vectors::mcols(gtf)), type == "gene" & gene_name == gene)]
  reads <- loadbamgene(bam = bam,
                       gene = S4Vectors::mcols(region)$gene_id,
                       region = region,
                       txdb = txdb,
                       exonByGene = EBG,
                       exonIntersect = exonIntersect,
                       junctionIntersect = junctionIntersect,
                       ignore.strand = FALSE,
                       isSupplementaryAlignment = isSupplementaryAlignment,
                       isSecondaryAlignment = isSecondaryAlignment,
                       isDuplicate = isDuplicate,
                       mapqFilter = mapqFilter,
                       sep = sep)
  reads <- reads[S4Vectors::mcols(reads)$CB %in% Cells]

  if(!is.null(groupBy)) {
    cov <- split(colnames(object), SummarizedExperiment::colData(object)[, groupBy])
    reads <- lapply(cov, function(x) {
      reads[S4Vectors::mcols(reads)$CB %in% as.character(x)]
    })
  } else {
    reads <- list("Depth" = reads)
  }
  reads <- lapply(reads, function(x) {
    if(S4Vectors::runValue(GenomicAlignments::strand(region)) == "+") {
      Pos <- data.table::data.table(pos = IRanges::start(IRanges::ranges(x)), CB = S4Vectors::mcols(x)$CB)
      cov <- GenomicAlignments::coverage(GenomicAlignments::first(x))[[unique(as.character(GenomeInfoDb::seqnames(region)))]]
    } else {
      Pos <- data.table::data.table(pos = IRanges::end(IRanges::ranges(x)), CB = S4Vectors::mcols(x)$CB)
      cov <- GenomicAlignments::coverage(GenomicAlignments::first(x))[[unique(as.character(GenomeInfoDb::seqnames(region)))]]
    }
    if(nrow(Pos) > 0) {
      cov <- data.table::data.table(pos = Pos[, min(pos)]:Pos[, max(pos)],
                                    D = as.numeric(cov[Pos[, min(pos)]:Pos[, max(pos)]]))
    } else {
      Pos <- data.table::data.table(pos = NA, CB = NA)
      cov <- data.table::data.table(pos = NA, D = NA)
    }
    return(list(cov, Pos))
  })
  Pos <- lapply(reads, function(x) x[[2]])
  cov <- lapply(reads, function(x) x[[1]])
  Pos <- data.table::data.table(groupBy = rep(names(Pos), mapply(nrow, Pos)), do.call(rbind, Pos))
  cov <- data.table::data.table(groupBy = rep(names(cov), mapply(nrow, cov)), do.call(rbind, cov))
  if(!is.null(groupLevels)) cov[, groupBy := factor(groupBy, levels = groupLevels)]
  if(!is.null(groupLevels)) Pos[, groupBy := factor(groupBy, levels = groupLevels)]
  Pos <- Pos[, .(y = seq_len(.N)), .(groupBy, pos)]
  # Tss <- lapply(Pos[, unique(groupBy)], function(x) {
  #   CliffSite2(x = Pos[groupBy == x, pos], SS = SSs, direction = "start", minreads = 50, max.degradation.length = NULL)
  # })

  if(is.null(xlimit)) {
    xlimit <- c(min(c(GenomicRanges::start(region), TSS)), max(c(GenomicRanges::end(region), TSS)))
  }
  gi <- biovizBase::crunch(txdb, region)
  gi <- gi[(GenomicRanges::start(gi) >= min(xlimit) & GenomicRanges::start(gi) <= max(xlimit)) | (GenomicRanges::end(gi) >= min(xlimit) & GenomicRanges::end(gi) <= max(xlimit))]
  GenomicRanges::start(gi[GenomicRanges::start(gi) < min(xlimit)]) <- min(xlimit)
  GenomicRanges::end(gi[GenomicRanges::end(gi) > max(xlimit)]) <- max(xlimit)
  gi <- GenomicRanges::GRangesList(split(gi, gi$tx_name))
  names(gi) <- plyr::mapvalues(names(gi), with(as.data.frame(gtf), transcript_id[type == "transcript"]), with(as.data.frame(gtf), transcript_name[type == "transcript"]), warn_missing = F)
  if(!is.null(transcripts)) {
    if(any(transcripts %in% names(gi))) {
      gi <- gi[names(gi) %in% transcripts]
    }
  }
  gi <- gi[mapply(length, gi) > 0]
  id <- mapply(function(x) paste(sort(c(start(x), end(x))), collapse = ","), gi)
  if(any(duplicated(id))) {
    id <- split(names(gi), id); names(id) <- NULL
    gi <- lapply(id, function(x) {
      gi[[x[1]]]
    })
    names(gi) <- mapply(function(x) paste(x, collapse = "; "), id)
  }
  gi <- GenomicRanges::GRangesList(gi)
  p0 <- ggbio::autoplot(gi)@ggplot + theme_void()

  cov <- cov[pos %between% ggplot2::layer_scales(p0)$x$range$range]
  Pos <- Pos[pos %between% ggplot2::layer_scales(p0)$x$range$range]
  cov <- cov[groupBy %in% cov[, max(D), groupBy][V1 >= MinDepth, groupBy]]
  Pos <- Pos[groupBy %in% Pos[, max(y), groupBy][V1 >= MinDepth, groupBy]]

  if(is.null(point.color)) {
    point.color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))(length(cov[, unique(groupBy)]))
  } else {
    if(length(point.color) == 1) {
      if(point.color %in% row.names(RColorBrewer::brewer.pal.info)) {
        point.color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info[point.color, "maxcolors"], point.color))(length(cov[, unique(groupBy)]))
      } else {
        e <- tryCatch(col2rgb(point.color), error = function(e) e)
        if(is(e, "simpleError")) {
          message("invalid color name, using grey to replace.")
          point.color <- rep("grey", length(cov[, unique(groupBy)]))
        } else {
          point.color <- rep(point.color, length(cov[, unique(groupBy)]))
        }
      }
    } else {
      if(length(point.color) != length(cov[, unique(groupBy)])) {
        message("invalid length of point.color, using RColorBrewer Set1 to replace.")
        point.color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))(length(cov[, unique(groupBy)]))
      }
    }
  }

  if(is.null(groupBy)) {
    if(is.null(Group.name.angle)) Group.name.angle <- 90
    ggplot(Pos) +
      geom_point(aes(x = pos, y = y, colour = groupBy), size = point.size, shape = point.shape) +
      scale_x_continuous(limits = ggplot2::layer_scales(p0)$x$range$range) +
      scale_y_continuous(n.breaks = 3) +
      scale_color_manual(values = point.color) +
      theme_light(base_size = base_size) +
      theme(legend.position = "none",
            axis.title = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.line.y = element_line(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.y.left = element_text(angle = Group.name.angle, colour = "black"),
            strip.placement = "outside",
            strip.background = element_rect(fill = NA)) -> p1
  } else {
    if(is.null(Group.name.angle)) Group.name.angle <- 0
    if(free_y) {
      ggplot(Pos) +
        geom_point(aes(x = pos, y = y, colour = groupBy), size = point.size, shape = point.shape) +
        facet_wrap(~ groupBy, ncol = 1, strip.position = "left", scales = "free_y") +
        scale_x_continuous(limits = ggplot2::layer_scales(p0)$x$range$range) +
        scale_y_continuous(n.breaks = 3) +
        scale_color_manual(values = point.color) +
        theme_light(base_size = base_size) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(),
              axis.line.y = element_line(),
              axis.ticks.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.y.left = element_text(angle = Group.name.angle, colour = "black"),
              strip.placement = "outside",
              strip.background = element_rect(fill = NA)) -> p1
    } else {
      ggplot(Pos) +
        geom_point(aes(x = pos, y = y, colour = groupBy), size = point.size, shape = point.shape) +
        facet_wrap(~ groupBy, ncol = 1, strip.position = "left") +
        scale_x_continuous(limits = ggplot2::layer_scales(p0)$x$range$range) +
        scale_y_continuous(n.breaks = 3) +
        scale_color_manual(values = point.color) +
        theme_light(base_size = base_size) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(),
              axis.line.y = element_line(),
              axis.ticks.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.y.left = element_text(angle = Group.name.angle, colour = "black"),
              strip.placement = "outside",
              strip.background = element_rect(fill = NA)) -> p1
    }
  }

  if(!is.null(TSS)) {
    p0 <- p0 + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
    # p1 <- p1 + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
  }

  if(is.null(rel_height)) {
    rh <- c(cov[, length(unique(groupBy))], max(0.5, length(gi)/10))
    rh <- rh / max(rh)
    rh[rh < 0.2] <- 0.2
    # ggbio::tracks(p1, p4@ggplot, heights = rh, padding = unit(0, "lines")) + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
    p1 / p0 + patchwork::plot_layout(ncol = 1, heights = rh)
  } else {
    # ggbio::tracks(p1, p4@ggplot, heights = c(1, rel_height), padding = unit(0, "lines")) + geom_vline(xintercept = TSS, colour = line.color, lty = line.type)
    p1 / p0 + patchwork::plot_layout(ncol = 1, heights = c(1, rel_height))
  }
}

