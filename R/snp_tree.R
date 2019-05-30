#' CART for Final SNP Selection
#'
#' Creates classification and regression trees for final SNP selection using
#' rpart. The function uses the dataset list from the \code{elasticnet}
#' function in this package as input and returns a list of trees and the
#' SNPs where the trees split. The trees are pruned using the 1-se rule.
#'
#' @param data_list Object obtained from the \code{elasticnet} function.
#' @param y A numeric vector of phenotypes for the data.
#'
#' @return \code{snp_tree} returns a list for each dataset in \code{data_list}.
#' Each of the sublists
#' contains the pruned tree created by rpart, a vector containing the SNPs where
#' the tree splits, and the alpha value
#' associated with the elastic net model used to create the dataset. Each
#' item is called using the list number
#' and name, such as results[[1]]$tree and results[[1]]$snps.
#'
#' \item{tree}{Tree produced by rpart. Tree is pruned using the 1-se rule.}
#' \item{snps}{Vector of SNPs where the tree splits. }
#' \item{alpha}{Alpha value from the elastic net model used to trim the
#' input data.}
#'
#' @export
#'
#'

snp_tree <- function(data_list, y) {

  if(!is.numeric(y)) stop("y vector must be numeric")
  if(length(unique(y)) == 2) y <- as.factor(y)

  tree <- vector("list", length(data_list))

  for(i in 1:length(data_list)) {

    if(dim(as.data.frame(data_list[[i]]$data))[2] > 1) {

      tmp_tree <- rpart::rpart(y ~ as.matrix(data_list[[i]]$data),
                               control = rpart::rpart.control(cp = 0.0, minsplit = 2))

      se1 <- tmp_tree$cptable[tmp_tree$cptable[, 4] == min(tmp_tree$cptable[, 4]), ]
      if(is.numeric(se1)) se1 <- data.frame(t(se1))
      se1.cp <- se1[1, 4] + se1[1, 5]
      cp <- tmp_tree$cptable[tmp_tree$cptable[, 4] < se1.cp, 1][1]

      snp_ <- as.matrix(data_list[[i]]$data)
      tree[[i]]$tree <- rpart::rpart(y ~ snp_, control = rpart::rpart.control(cp = cp))

      tst <- as.character(tree[[i]]$tree$frame$var[!grepl("leaf", tree[[i]]$tree$frame$var)])
      tree[[i]]$snps <- sub(".*)", "", tst)

      tree[[i]]$alpha <- data_list[[i]]$alpha

    } else {

      tree[[i]]$alpha <- data_list[[i]]$alpha

    }
  }

  tree
}
