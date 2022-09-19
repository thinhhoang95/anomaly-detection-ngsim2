oldPar <- par(no.readonly = TRUE)
set.seed(1)
### simulate data (one-dimensional domains)
sim <- simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
                       M = 5, eFunType = "Poly", eValType = "linear", N = 100)
plot(sim$simData)
# MFPCA based on univariate FPCA
uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"), list(type = "uFPCA")))
summary(uFPCA)
plot(uFPCA) # plot the eigenfunctions as perturbations of the mean
scoreplot(uFPCA) # plot the scores
# MFPCA based on univariate spline expansions
splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                                          list(type = "splines1D", k = 10)), fit = TRUE) # calculate reconstruction, too
summary(splines)
plot(splines) # plot the eigenfunctions as perturbations of the mean
scoreplot(splines) # plot the scores
