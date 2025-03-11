#' scATSDataSet

setClass("scATSDataSet",
         contains = "RangedSummarizedExperiment",
         slots = representation(
           reductions = "list",
           project.name = "ANY"
         ))

scATSDataSet <- function(se, reductions = NULL, project.name = NULL, ...) {
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }
  # if (!identical(SummarizedExperiment::assayNames(se), c("counts", "expected.counts", "psi"))) {
  #   stop("counts, expected.counts, and psi slots are necessary")
  # }
  if (!all(is.element(c("counts", "psi"), SummarizedExperiment::assayNames(se)))) {
    stop("counts, and psi slots are necessary")
  }

  object <- new("scATSDataSet", se, reductions = reductions, project.name = project.name, ...)
  # stash the package version
  S4Vectors::metadata(object)[["version"]] <- packageVersion("scATS")
  # validObject(object)
  return(object)
}

#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment assay
#' @export
setGeneric("expected.counts", function(object, cells = NULL, TSSs = NULL) standardGeneric("expected.counts"))

setMethod(f = "expected.counts",
          signature = "scATSDataSet",
          definition = function(object, cells = NULL, TSSs = NULL) {
            if(!is.null(cells)) stopifnot(all(is.element(cells, colnames(object))))
            if(!is.null(TSSs)) stopifnot(all(is.element(TSSs, row.names(object))))
            cells <- cells %||% colnames(object)
            TSSs <- TSSs %||% row.names(object)
            SummarizedExperiment::assay(object[TSSs, cells], "expected.counts")
          })

#' @export
setGeneric("psi", function(object, cells = NULL, TSSs = NULL, groupBy = NULL) standardGeneric("psi"))

setMethod(f = "psi",
          signature = "scATSDataSet",
          definition = function(object, cells = NULL, TSSs = NULL, groupBy = NULL) {
            if(!is.null(cells)) stopifnot(all(is.element(cells, colnames(object))))
            if(!is.null(TSSs)) stopifnot(all(is.element(TSSs, row.names(object))))
            cells <- cells %||% colnames(object)
            TSSs <- TSSs %||% row.names(object)
            if(!is.null(groupBy)) {
              stopifnot(length(groupBy) == 1)
              stopifnot(groupBy %in% colnames(SummarizedExperiment::colData(object)))
            }
            res <- SummarizedExperiment::assay(object[TSSs, cells], "psi")
            if(!is.null(groupBy)) {
              res <- merge(data.table::as.data.table(t(res), keep.rownames = "Cells"), data.table::data.table(Cells = row.names(SummarizedExperiment::colData(object)), groupBy = SummarizedExperiment::colData(object)[, groupBy]), by = "Cells")
              res <- data.table::melt.data.table(res, id.vars = c("Cells", "groupBy"), variable.name = "TSS", value.name = "PSI", variable.factor = FALSE)
              res <- res[, .(Cells = .N, N = base::sum(!is.na(PSI)),
                             n = base::sum(PSI > 0, na.rm = T),
                             mean = base::mean(PSI, na.rm = T),
                             sd = stats::sd(PSI, na.rm = T),
                             se = stats::sd(PSI, na.rm = T)/base::sqrt(base::sum(!is.na(PSI))),
                             ci = stats::qt(0.95/2 + 0.5, base::sum(!is.na(PSI)) - 1) * stats::sd(PSI, na.rm = T)/base::sqrt(base::sum(!is.na(PSI))),
                             median = stats::median(PSI, na.rm = T),
                             Q1 = stats::quantile(PSI, probs = 0.25, na.rm = T),
                             Q3 = stats::quantile(PSI, probs = 0.75, na.rm = T),
                             mad = stats::mad(PSI, na.rm = T),
                             iqr = stats::IQR(PSI, na.rm = TRUE)
              ), .(TSS, groupBy)]

              genes <- mapply(function(x) x[1], strsplit(TSSs, "@"))

              if(identical(TSSs, row.names(object))) {
                TSSss <- row.names(object)
              } else {
                TSSss <- row.names(object)[unique(do.call(c, lapply(FUN = function(x) grep(x, row.names(object)), genes)))]
              }

              mat_p <- SummarizedExperiment::assay(object[TSSss, cells], "psi")
              mat_n <- SummarizedExperiment::assay(object[TSSss, cells], "counts")
              mat_n <- reshape2::melt(mat_n, value.name = "n")
              mat_p <- reshape2::melt(mat_p, value.name = "PSI")
              if(identical(mat_n[, 1:2], mat_p[, 1:2])) {
                mat <- cbind(mat_n, PSI = mat_p[, 3])
              } else {
                mat <- merge(mat_n, mat_p, by = c("Var1", "Var2"))
              }
              mat <- data.table::as.data.table(mat)

              mat$groupBy <- SummarizedExperiment::colData(object[TSSss, cells])[as.character(mat$Var2), groupBy]
              mat$gene_id <- SummarizedExperiment::rowData(object[TSSss, cells])[as.character(mat$Var1), "gene_id"]
              data.table::setkey(mat, gene_id, Var2)
              N <- suppressWarnings(mat[, .(N = max(n/PSI, na.rm = T)), .(gene_id, Var2)])
              N[is.infinite(N), N := 0]
              mat <- data.table::merge.data.table(mat, N, by = c("gene_id", "Var2"))
              mat <- mat[, .(PseudobulkPSI = sum(n, na.rm = TRUE) / sum(N, na.rm = TRUE)), .(Var1, groupBy)]
              res <- as.data.frame(merge(res, mat, by.x = c("TSS", "groupBy"), by.y = c("Var1", "groupBy")))
            }
            return(res)
          })

#' @export
setGeneric("counts", function(object, cells = NULL, TSSs = NULL) standardGeneric("counts"))

setMethod(f = "counts",
          signature = "scATSDataSet",
          definition = function(object, cells = NULL, TSSs = NULL) {
            if(!is.null(cells)) stopifnot(all(is.element(cells, colnames(object))))
            if(!is.null(TSSs)) stopifnot(all(is.element(TSSs, row.names(object))))
            cells <- cells %||% colnames(object)
            TSSs <- TSSs %||% row.names(object)
            SummarizedExperiment::assay(object[TSSs, cells], "counts")
          })

#'

setClass("SimpleSeurat",
         contains = "list",
         slots = representation(
           Cell = "character",
           Gene = "character",
           meta.data = "ANY",
           reductions = "list"
         ))

setMethod(
  f = 'show',
  signature = 'SimpleSeurat',
  definition = function(object) {
    cat(" No. genes:", length(object@Gene), "\n", "No. cells: ", length(object@Cell), "\n")
  }
)

#' @export
#'
SimplifySeurat <- function(object) {
  new("SimpleSeurat", Cell = colnames(object), Gene = row.names(object), meta.data = object@meta.data, reductions = object@reductions)
}

