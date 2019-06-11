#' Subset Dataset Based on Filter Results
#'
#' \code{filter_subset} uses the original genotypes and the vector returned by
#' \code{snp_filter} to create a new dataset that contains only the variables
#' that pass through the user-specified threshold for p-values or distance
#' correlations. If a window is used,  the returned dataset will
#' include all of the SNPs in a window for a SNP that meets the threshold.
#'
#' @param x A matrix or dataset containing the genotypes for each
#' observation. The number of variables must match the number of
#' elements in the \code{results} vector.
#' @param results A vector of p-values or distance correlations,
#' preferably the output vector from the \code{snp_filter} function.
#' The length of \code{results} must be equal to the number of variables
#' in \code{x}.
#' @param type The type of measure in the results vector. The
#' default is "dc" for distance correlations. If a different
#' string is used it is assumed that the vector contains p-values.#'
#' @param threshold This value can either be decimal less than 1 or
#' a positive integer. If a number less than 1 and if
#' the results are p-values, all SNPs with a p-value less than
#' the threshold will be retained. If the results are distance
#' correlations, all SNPs with a value greater than the threshold
#' will be retained. A positive integer indicates how many SNP
#' should pass through the filter.
#' @param window This is the size of the window used in the initial filter.
#' If a window was used in the \code{snp_filter} function. The dataset that
#' is returned will contain the selected SNP and the SNPs in the
#' window around the SNP.
#'
#' @return Returns a dataset that contains the SNPs
#' that pass through the user specified threshold. The output can
#' be used directly as an input into the \code{elasticnet} function in this
#' package.
#'
#' @export
#'

filter_subset <- function(x, results, type = "dc", threshold, window = 1) {

  if(ncol(x) != length(results)) stop("Number of results is not equal
                                        to the number of variables in x.")
  if(window < 1) stop("Invalid value of window. Use a positive integer")
  if(window > ncol(x)) stop("Invalid value of window. Use a positive integer
                            less than the number of SNPs in x.")

  if(threshold > 1) {
    if(type == "dc") {
      threshold <- results[order(results, decreasing = TRUE)][threshold]
    } else {
      threshold <- results[order(results)][threshold]
    }
  }

  if(type == "dc") {
    results <- -1 * results
    threshold <- -1 * threshold
  }

  if(window == 1) {

    dat <- x[, results <= threshold]

  } else {

    window <- round(window / 2)
    keep <- NULL
    keep2 <- 1:ncol(x)
    keep2 <- keep2[results <= threshold]

    for(i in 1:length(keep2)) {
      rng <- (keep2[i] - window):(keep2[i] + window)
      rng <- rng[rng > 0]
      rng <- rng[rng <= ncol(x)]
      keep <- c(keep, rng)
    }

    keep <- unique(keep)
    keep <- keep[order(keep)]
    dat <- x[, keep]
  }

  dat
}
