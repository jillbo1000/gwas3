#' First Phase Filtering
#'
#' \code{snp_filter} performs the preliminary filter calculations. Filters can
#' be either distance correlation or linear/logistic regression. The function
#' will determine if the phenotype is binary and will select the correct
#' model type for the user when the linear/logistic regression option
#' is selected. Although the documentation is written in terms of GWAS
#' data, the function can be used on any dataset with numeric variables.
#'
#' @param x A matrix or dataset containing the genotypes for each observation.
#' @param y A numeric vector containing the phenotypes. Can be binary or continuous.
#' @param type The type of filter to perform. The default is distance
#' correlation. If any string is insterted other than "dc", either linear
#' regression or logistic regression will be used as the filter.
#' @param window The size of the window, or neighborhood, around each SNP that
#' will be used to compute the linear or logistic regression model. If a window of
#' 10 is chosen, the five SNPs on either side of the selected SNPs will be used
#' to compute the distance correlation or linear model. The default is single
#' SNP analysis.
#' @param fdr If TRUE, the false-discovery rate correction for the p-value is
#' returned instead of the raw p-values. If type = "dc", this is ignored.
#'
#' @return Returns a vector that contains the distance correlations
#' or p-values for each observation. If window = 1, the p-value is the
#' p-value for the coefficient of the SNP. If window > 1 the p-value is
#' from the likelihood ratio test for the model.
#'
#' @export
#'
snp_filter <- function(x, y, type = "dc", window = 1, fdr = TRUE) {

  if(!is.numeric(y)) stop("y vector must be numeric")
  if(window < 1) stop("Invalid value of window. Use a positive integer")
  if(window > ncol(x)) stop("Invalid value of window. Use a positive integer
                            less than the number of SNPs in x.")

  results <- rep(NA, ncol(x))
  window <- round(window / 2)

  if(type == "dc") {
    for(i in 1:length(results)) {
      rng <- (i - window):(i + window)
      rng <- rng[rng > 0]
      rng <- rng[rng <= ncol(x)]
      tmp <- x[, rng]
      try(results[i] <- energy::dcor(tmp, y))
    }
  } else {

    if(length(unique(y)) == 2) {
      y <- as.factor(y)
      for(i in 1:length(results)) {
        rng <- (i - window):(i + window)
        rng <- rng[rng > 0]
        rng <- rng[rng <= ncol(x)]
        tmp <- x[, rng]
        fit <- stats::glm(y ~ as.matrix(tmp), family = "binomial")
        if(window < 1) {
          try(results[i] <- stats::coef(summary(fit))[2,4])
        } else {
          try(results[i] <- lmtest::lrtest(fit)$Pr[2])
        }
      }
    } else {

      for(i in 1:length(results)) {
        rng <- (i - window):(i + window)
        rng <- rng[rng > 0]
        rng <- rng[rng <= ncol(x)]
        tmp <- x[, rng]
        fit <- stats::glm(y ~ as.matrix(tmp), family = "gaussian")
        if(window < 1) {
          try(results[i] <- stats::coef(summary(fit))[2,4])
        } else {
          try(results[i] <- lmtest::lrtest(fit)$Pr[2])
        }
      }
    }
  }

  if(fdr == TRUE & type != "dc") {results <- stats::p.adjust(results, method = "BH")}
  results
}
