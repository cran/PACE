FPCA = function(y, tt, p = setOptions()){

   X = PCA(y,tt, bwmu = p$bwmu, bwmu_gcv = p$bwmu_gcv, bwxcov = p$bwxcov ,bwxcov_gcv = p$bwxcov_gcv,
               ntest1 = p$ntest1, ngrid1 = p$ngrid1, selection_k = p$selection_k, FVE_threshold = p$FVE_threshold,
               maxk = p$maxk, control = p$control, regular = p$regular, error = p$error, ngrid = p$ngrid,
               method = p$method, shrink = p$shrink, newdata = p$newdata, kernel = p$kernel,
               numBins = p$numBins, yname = p$yname, screePlot = p$screePlot, designPlot = p$designPlot, rho = p$rho,
               verbose = p$verbose)
   X
}