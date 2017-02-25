Conditional Analysis of Multiple Traits Based on Marginal GWAS Summary Statistics

#Joint model for L SNPs and K traits. (JointSum)
#Model 1: trait 1 ~ L SNPs + (K-1) traits.
#Model 2: trait 1 ~ L SNPs.

#INPUT
#B1: marginal effects on trait 1. It should be a L*1 matrix.
#S1: standard errors for marginal effects on trait 1. L*1.
#B2: marginal effects on traits to adjust for. It should be a L*(K-1) matrix.
#S2: standard errors that correspond to B2. L*(K-1).
#N: sample sizes for each coefficients in B1, B2. It should be a L*K matrix where the first column corresponds to the sample sizes for marginal effects on trait 1.
#XX: estimated covariance matrix for SNPs. Scaling this matrix will not affect the results.
#YY0: estimated correlation matrix for traits.
#adj_Y: if it is 0, adjust for SNPs only. Otherwise adjust for both SNPs and traits.
#lam: a modifying parameter in [0,1). It is used only if adj_Y=1.

#OUTPUT
#beta: coefficient estimates.
#cov: covariance matrix for coefficients.
#pvalue: p-values for coefficients.
#sigma2: estimated MSE.
