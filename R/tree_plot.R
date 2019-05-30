#' Plot Trees from Final SNP Selection
#'
#' Plots the trees obtained using \code{snp_tree}.
#'
#' @importFrom graphics "par"
#'
#' @param tree Object obtained from the \code{snp_tree} function.
#' @param nrow,ncol Integers that specify the number of rows and
#' columns in the plotting screen.
#' @param palette User specified palette using options from \code{rpart.plot}.
#'
#' @return \code{plot_tree} plots the trees obtained from \code{snp_tree} using
#' the \code{rpart.plot} package. If a plot is missing for a value of
#' alpha, it is because one or no variables were retained from the elastic net
#' model.
#'
#' @export
#'
#'

tree_plot <- function(tree, nrow = 2, ncol = 4, palette = "BlGnYl") {
  oldpar <- graphics::par()

  par(mfrow = c(nrow, ncol))

  for(i in 1:length(tree)) {

    if(!is.null(tree[[i]]$tree)) {
      if(tree[[i]]$tree$method == "class") {
        rpart.plot::rpart.plot(tree[[i]]$tree, extra = 102, under = TRUE,
                               box.palette = palette,
                               main = paste("alpha = ", tree[[i]]$alpha, sep = ""))
      } else {
        rpart.plot::rpart.plot(tree[[i]]$tree, under = TRUE, box.palette = palette,
                               main = paste("alpha = ", tree[[i]]$alpha, sep = ""))
      }
    }
  }

  graphics::par(oldpar)
}
