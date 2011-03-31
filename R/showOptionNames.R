showOptionNames = function(){

   cat("Here is a list of optional input argument names for the PCA():\n")
   op_names = c("bwmu", "bwmu_gcv", "bwxcov", "bwxcov_gcv", "ntest1", "ngrid1",
                "selection_k", "FVE_threshold", "maxk", "control", "regular", 
                "error", "ngrid", "method", "shrink", "newdata", "kernel",
                "numBins", "yname", "screePlot", "designPlot", "rho", "verbose")
   print(op_names)
}

