setOptions = function(bwmu = 0, bwmu_gcv = 1, bwxcov = c(0,0), bwxcov_gcv = 1, 
                      ntest1 = 30, ngrid1 = 30, selection_k = "BIC1", FVE_threshold = 0.85,
                      maxk = 20, control = "auto", regular = 0, error = 1, ngrid = 51,
                      method = "CE", shrink = 0, newdata = NULL, kernel = "gauss", 
                      numBins = NULL, yname = NULL, screePlot = 0, designPlot = 0, rho = "cv", verbose = "on")
{

  list(bwmu = bwmu, bwmu_gcv = bwmu_gcv, bwxcov = bwxcov, bwxcov_gcv = bwxcov_gcv,
       ntest1 = ntest1, ngrid1 = ngrid1, selection_k = selection_k, FVE_threshold = FVE_threshold,
       maxk = maxk, control = control, regular = regular, error = error, ngrid = ngrid, 
       method = method, shrink = shrink, newdata = newdata, kernel = kernel, 
       numBins = numBins, yname = yname, screePlot = screePlot, designPlot = designPlot, rho = rho, 
       verbose = verbose)

}
