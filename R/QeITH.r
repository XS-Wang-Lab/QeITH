#' quantifying intratumor immune heterogeneity
#'
#' @param data A data.frame of mRNA expression profile, cell types annotation data, or the proportions of cell types across samples.
#' @param type A character vector indicating the type of "bulk", "single-cell",'spatial'.
#'
#' @returns A data frame containing sample IDs and intratumor heterogeneity scores -- ITH score
#' @export
#'
#' @examples
#' # Example usage:
#' # result_bulk <- QeITH(data, type = bulk")
#' # result_sc <- QeITH(data, type = "single-cell")
#' # result_sc <- QeITH(data, type = "spatial")

QeITH <- function(data, type = c('bulk', 'single-cell','spatial')) {

  # Parameter validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input 'data' must be a matrix or data.frame")
  }

  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("Input data is empty or has zero dimensions")
  }

  if (missing(type)) {
    stop("The type parameter must be specified. You can choose: 'bulk', 'single-cell', or 'spatial'")
  }

  type <- match.arg(type)

  # Process based on input type
  if (type == "bulk") {

    load(system.file("data","cellular_states.RData", package = "QeITH"))

    # Check if GSVA package is available
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      stop("Package 'GSVA' is required but not installed.")
    }

    geneset_list <- lapply(marker, function(x) x[!is.na(x)])
    rna_seq = as.matrix(data)

    # Check gene overlap between data and gene sets
    common_genes <- intersect(rownames(rna_seq), unlist(geneset_list))
    if (length(common_genes) < 10) {
      warning(sprintf("Only %d genes overlap between data and gene sets", length(common_genes)))
    }

    results <- GSVA::gsva(rna_seq, geneset_list,
                          method="ssgsea",
                          ssgsea.norm=TRUE,
                          verbose=FALSE)


    results = t(results)
    results[results < 0] <- 0

    # Calculate row sums for normalization
    row_sums <- rowSums(results, na.rm = TRUE)
    row_sums[row_sums == 0] <- 1  # 或者设为1
    ratio_matrix <- sweep(results, 1, row_sums, FUN = "/")

  }  else if (type == "single-cell") {
    # Validate required columns for single-cell mode
    if (!all(c("Samples", "Celltype") %in% colnames(data)))     {
      stop("For single-cell mode, data must contain 'Samples' and 'Celltype' columns")
    }

    # Check required packages for single-cell mode
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required for single-cell mode.")
    }
    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("Package 'tidyr' is required for single-cell mode.")
    }

    library(dplyr, quietly = TRUE)
    # Count cells per sample and cell type
    result <- data %>%
      dplyr::group_by(Samples, Celltype) %>%
      dplyr::summarise(count = n(), .groups = 'drop')

    # Convert to wide format (sample × cell type matrix)
    result2 <- result %>%
      tidyr::pivot_wider(names_from = Celltype, values_from = count, values_fill = list(count = 0)) %>%
      as.data.frame()

    # Convert to proportion matrix
    ratio_matrix <- result2
    rownames(result2) <- result2$Samples
    result2 <- result2[, -1, drop = FALSE]

    cell_counts <- as.matrix(result2)
    row_sums <- rowSums(cell_counts, na.rm = TRUE)
    ratio_matrix <- cell_counts / row_sums

  } else if (type == "spatial") {
    # For spatial transcriptomics data, use the cell proportion matrix directly
    ratio_matrix <- as.matrix(data)
  }

  # Check the range of cell proportions
  if (any(ratio_matrix < 0, na.rm = TRUE) || any(ratio_matrix > 1, na.rm = TRUE)) {
    # 给出一个提示信息，说明正在调整
    warning(
      "Some values were slightly outside [0,1] range due to numerical precision. ",
      "They have been clipped to the nearest boundary (0 or 1)."
    )
    # 将小于0的值设置为0，大于1的值设置为1
    ratio_matrix[ratio_matrix < 0] <- 0
    ratio_matrix[ratio_matrix > 1] <- 1
  }

  # Normalize rows to ensure they sum to 1
  row_sums <- rowSums(ratio_matrix, na.rm = TRUE)
  needs_norm <- abs(row_sums - 1) > 1e-6

  if (any(needs_norm)) {
    warning(sprintf("%d rows normalized to sum to 1", sum(needs_norm)))
    ratio_matrix[needs_norm, ] <- ratio_matrix[needs_norm, ] / row_sums[needs_norm]
  }

  # Calculate ITH score
  calculate_entropy <- function(x) {
    x <- x[x > 0 & !is.na(x)]  # Remove zeros
    if (length(x) == 0) return(NA)
    x <- x / sum(x)
    -sum(x * log2(x))
  }

  scores <- apply(ratio_matrix, 1, calculate_entropy)

  # Output results of ITH scores for all samples
  result <- data.frame(
    Sample = rownames(ratio_matrix),
    ITH_Score = as.numeric(scores),
    stringsAsFactors = FALSE
  )

  return(result)
}

