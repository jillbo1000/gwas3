library(testthat)
library(gwas3)
library(scales)

test_check("gwas3")

mouse <- read.csv("C:/Users/jflun/Dropbox/USU/Advanced Regression/Project/Mouse Data/Mouse.csv")
chr <- read.csv("C:/Users/jflun/Dropbox/USU/Advanced Regression/Project/Mouse Data/Mouse Position.csv")

dim(chr)
head(chr)
tail(chr)

dim(mouse)
head(mouse[, 1:10])

m_geno <- mouse[, -(1:2)]
# m_geno <- m_geno[, 1:5000] # Reduce size of mouse for faster computation
head(m_geno[, 1:10])
write.csv(m_geno, "C:/Users/jflun/Dropbox/Dissertation/LASSO_EN_Project/Data/mouse/mouse_geno.csv", row.names = FALSE)

m_sim <- gwas_sim(m_geno, 10, 0.5)
m_sim

summary(m_sim$phenotype)
m_sim$functional_snps

# Plot for dc
mgeno_dc <- snp_filter(m_geno, m_sim$phenotype, type = "dc", window = 1)
plot(1:length(mgeno_dc), mgeno_dc, pch = 16, col = alpha("gray", 0.5), cex = 0.7)
# points(m_sim$functional_snps, rep(0.3, length(m_sim$functional_snps)), col = "red", pch = 16, cex = 0.7)

effect_dc <- abs(m_sim$effect) / max(abs(m_sim$effect)) * max(mgeno_dc)
lines(m_sim$functional_snps, effect_dc, col = "blue")

for(i in 1:length(m_sim$functional_snps)) {
  # abline(v = m_sim$functional_snps[i], col = "red")
  lines(c(m_sim$functional_snps[i], m_sim$functional_snps[i]), c(0, effect_dc[i]), col = "blue")
}

# plot for linear regression
mgeno_lm <- snp_filter(m_geno, m_sim$phenotype, type = "lm", window = 1)
plot(1:length(mgeno_lm), -log10(mgeno_lm), pch = 16, col = alpha("gray", 0.5), cex = 0.7)
# points(m_sim$functional_snps, rep(0.3, length(m_sim$functional_snps)), col = "red", pch = 16, cex = 0.7)

effect_lm <- abs(m_sim$effect) / max(abs(m_sim$effect)) * max(-log10(mgeno_lm))
lines(m_sim$functional_snps, effect_lm, col = "blue")

for(i in 1:length(m_sim$functional_snps)) {
  # abline(v = m_sim$functional_snps[i], col = "red")
  lines(c(m_sim$functional_snps[i], m_sim$functional_snps[i]), c(0, effect_lm[i]), col = "blue")
}


# plot for linear regression
mgeno_lm_10 <- snp_filter(m_geno, m_sim$phenotype, type = "lm", window = 10)
plot(1:length(mgeno_lm_10), -log10(mgeno_lm_10), pch = 16, col = alpha("gray", 0.5), cex = 0.7)
# points(m_sim$functional_snps, rep(0.3, length(m_sim$functional_snps)), col = "red", pch = 16, cex = 0.7)

effect_lm <- abs(m_sim$effect) / max(abs(m_sim$effect)) * max(-log10(mgeno_lm_10))
lines(m_sim$functional_snps, effect_lm_10, col = "blue")

for(i in 1:length(m_sim$functional_snps)) {
  # abline(v = m_sim$functional_snps[i], col = "red")
  lines(c(m_sim$functional_snps[i], m_sim$functional_snps[i]), c(0, effect_lm[i]), col = "blue")
}


plot(1:length(mgeno_lm_10), mgeno_lm_10, pch = 16, col = alpha("gray", 0.5), cex = 0.7)


par(mfrow = c(1, 3))
plot(1:length(mgeno_lm), -log10(mgeno_lm), pch = 16, col = alpha("black", 0.5), cex = 1)
effect_lm <- abs(m_sim$effect) / max(abs(m_sim$effect)) * max(-log10(mgeno_lm))
for(i in 1:length(m_sim$functional_snps)) {
  lines(c(m_sim$functional_snps[i], m_sim$functional_snps[i]), c(0, effect_lm[i]), col = "red")
}

plot(1:length(mgeno_lm), -log10(mgeno_lm), pch = 16, col = alpha("black", 0.5), cex = 1)
effect_lm <- abs(m_sim$effect) / max(abs(m_sim$effect)) * max(-log10(mgeno_lm))
for(i in 1:length(m_sim$functional_snps)) {
  lines(c(m_sim$functional_snps[i], m_sim$functional_snps[i]), c(0, effect_lm[i]), col = "red")
}

plot(1:length(mgeno_lm_10), -log10(mgeno_lm_10), pch = 16, col = alpha("black", 0.5), cex = 1)
effect_lm_10 <- abs(m_sim$effect) / max(abs(m_sim$effect)) * max(-log10(mgeno_lm_10))
for(i in 1:length(m_sim$functional_snps)) {
  lines(c(m_sim$functional_snps[i], m_sim$functional_snps[i]), c(0, effect_lm_10[i]), col = "red")
}

