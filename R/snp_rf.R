#' Random Forests for Final SNP Selection
#'
#' Computes the importance of each SNP in the input datasets and creates
#' variable importance plots for final SNP selection using
#' random forests. The function
#' uses the dataset list returned by the \code{elasticnet} function
#' as input.
#'
#' @param data_list Object obtained from the \code{elasticnet} function in this
#' package.
#' @param y A numeric vector of phenotypes for the data. Can be either
#' continuous or binary.
#'
#' @return \code{snp_rf} returns a list for each dataset in
#' \code{data_list}. Each element of the list contains
#' the variable importance measures, variable importance plot,
#' and the alpha value associated with the \code{elasticnet} model used
#' to create the dataset. Each item is called using the list number
#' and name, such as results[[1]]$vi
#'
#' \item{vi}{Vector that includes the SNPs and their importance.}
#' \item{vi_plot}{Variable importance plot.}
#' \item{alpha}{Alpha value from the elastic net model used to trim the
#' input data.}
#'
#' @export
#'

snp_rf <- function(data_list, y) {

  if(!is.numeric(y)) stop("y vector must be numeric")
  if(length(unique(y)) == 2) y <- as.factor(y)

  rf <- vector("list", length(data_list))

  # The following variable assignment is made only to appease RMD CHECK
  SNP <- Importance <- NULL

  for(i in 1:length(data_list)) {

    if(dim(as.data.frame(data_list[[i]]$data))[2] > 1) {

      tmp_rf <- randomForest::randomForest(data_list[[i]]$data, y)
      rf[[i]]$vi <- randomForest::importance(tmp_rf)

      tmp_best <- data.frame(SNP = row.names(randomForest::importance(tmp_rf)))
      tmp_best$Importance <- randomForest::importance(tmp_rf)
      limit <- min(c(15, nrow(tmp_best)))
      best15 <- tmp_best[order(tmp_best$Importance, decreasing = TRUE), ][1:limit, ]
      best15$SNP <- factor(best15$SNP, levels = rev(best15$SNP))
      rng <- (range(best15$Importance)[2] - range(best15$Importance)[1]) * 0.02

      rf[[i]]$vi_plot <- ggplot2::ggplot(best15, ggplot2::aes(x = SNP, y = Importance)) +
        ggplot2::geom_dotplot(binaxis = "y", method = "dotdensity",
                     binwidth = 1, dotsize = rng) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste("alpha = ", data_list[[i]]$alpha, sep = "")) +
        ggplot2::coord_flip()

      rf[[i]]$alpha <- data_list[[i]]$alpha

    } else {

      rf[[i]]$alpha <- data_list[[i]]$alpha
    }
  }

  rf
}
