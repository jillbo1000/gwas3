#' Second Phase of Filtering
#'
#' \code{elasticnet} uses the genotypes and phenotypes that pass through
#' the first phase of SNP filtering and computes LASSO and elastic net
#' models to further filter the SNPs. Lambda is determined using
#' cross-validation with the LASSO model. That value of lambda
#' is used for all elastic net
#' and LASSO models. A list containing datasets of genotypes are
#' returned for alpha values of 0, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75,
#' 0.9, and 1.0. Note that the model with alpha = 0 will retain all of
#' the SNPs passed into the function. This option is included so the user
#' can easily compare the performance of the final SNP selction methods without
#' the second stage filtering with models that have the additional filtering.
#' The final datasets contain the SNPs that have non-zero coefficients.
#'
#' @param x A matrix or dataframe containing the genotypes that passed
#' through the first phase of filtering.
#' @param y A numeric vector containing the phenotypes. The function will
#' determine if the response is binary or continuous.
#'
#' @return Returns a list where each member of the list is a dataset
#' that contains the genotypes for the SNPs that have non-zero
#' coefficients and the alpha value used in the elastic net model.
#' The list elements are numbered rather than named to make it easier
#' to use as inputs into other functions. For example, result[[1]]$data
#' contains the dataset and result[[1]]$alpha contains the value of alpha
#' used to compute the model.
#'
#' @export

elasticnet <- function(x, y) {

  if(length(unique(y)) == 2) {
    fam <- "binomial"
    y <- as.factor(y)
  } else {
    fam <- "gaussian"
    y <- as.numeric(as.character(y))
  }

  test <- glmnet::cv.glmnet(as.matrix(x), y, family = fam)
  alpha <- c(0, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
  dat <- vector("list", length(alpha))

  for(i in 1:length(alpha)) {
    test2 <- glmnet::glmnet(as.matrix(x), y, family = fam,
                            lambda = test$lambda.1se, alpha = alpha[i])
    b <- stats::coef(test2)[-1]
    x_en <- x[, abs(b) > 0]
    dat[[i]]$data <- x_en
    dat[[i]]$alpha <- alpha[i]

  }

  dat
}
