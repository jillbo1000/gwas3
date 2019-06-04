#' Simulate GWAS data with known functional SNPs
#'
#' \code{gwas_sim} takes a real set of genotypes and simulates a phenotype
#' with a known set of important SNPs. The user can select the number of
#' SNPs that will effect the phenotype and the heretability. The function
#' returns the phenotype, the list of SNPs that are important, a dataset
#' that contains the genotypes with the important SNPs removed, and some
#' other features that may be of interest. Note that the genotypes are
#' not changed. The phenotype is created to be impacted by specified SNPs.
#' This method preserves all of the physical characteristics of the
#' organisms genetics.
#'
#' Phenotypes are computed by randomly selecting an set of SNPs that
#' will be functional and creating random normal variables
#' with mean 0 and variance 1 to be
#' effects of the functional SNPs. The effects for
#' the rest of the SNPs are 0. The phenotype is y = effect * SNP + noise.
#' The noise is a random normal variables with mean = 0 and
#' variance = (1 - heretability) * var(effect * SNP) / heretability.
#'
#'
#' @param genotype A matrix containing the genotypes for an organism.
#' The rows represent the subjects and the columns are the SNPs.
#' @param num_snps The number of snps that will influence the phenotype.
#' @param heretability A number between 0 and 1 that is the proportion
#' variance in the phenotype explained by the genotype (R-squared).
#'
#' @return Returns a phenotype that is directly impacted by the specified
#' SNPs and other metrics that are useful in using and interpreting the
#' simulated data.
#'
#' \item{phenotype}{A numeric vector with the phenotype for the
#' simulated data. This is continuous, but a binary phenotype can
#' be created by the user from the returned phenotypes.}
#' \item{functional_snps}{Numeric vector containing the column number
#' of the SNPs that are impacting the phenotype.}
#' \item{geno_snp_removed}{A dataset of genotypes with the functional
#' SNPs removed.}
#' \item{estimated_heretability}{A measure of the heretability in the
#' simulated dataset. This should be close to the heretability that was
#' entered as an argument.}
#' \item{effect}{Effect measures used as a coefficient for each important
#' SNP.}
#'
#' @export
#'

gwas_sim <- function(genotype, num_snps = 10, heretability = 0.3) {

  snps <- sample(1:ncol(genotype), num_snps)
  snps <- snps[order(snps)]
  effect <- stats::rnorm(num_snps, 0, 1)
  ph <- rep(NA, nrow(genotype))

  for(i in 1:nrow(genotype)) {
    ph[i] <- sum(genotype[i, snps] * effect)
  }

  pheno_var <- (1 - heretability) * stats::var(ph) / heretability
  phenotype <- ph + stats::rnorm(nrow(genotype), 0, sqrt(pheno_var))

  sim <- NULL

  sim$phenotype <- phenotype
  sim$functional_snps <- snps
  sim$geno_snp_removed <- genotype[, -snps]
  sim$estimated_heretability <- stats::var(ph) / stats::var(phenotype)
  sim$effect <- effect

  sim

}

