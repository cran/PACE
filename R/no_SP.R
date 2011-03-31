no_SP = function(FVE, ngrid, no_optCopy, FVE_threshold = 0.85, yname = "y", verbose = "on"){

       fve2plot = FVE[1:ceiling(ngrid/2)]*100;
       fve2plot = fve2plot[!is.na(fve2plot)]
       x11()
       x = 0:length(fve2plot)
       y = c(0,fve2plot)
       plot(x, y, type = "p", lwd = 3, xlim = c(0, length(FVE)+1), ylim = c(0,101),
            main = paste("Fraction of variance explained by No. of PC \n(threshold = ", FVE_threshold, 
            ") for function ", yname, sep = ""), xlab = "No. of Principal Components", 
            ylab = "FVE (%)", font.lab = 2)
       
       lines(x,y, lty = 2, col = "red", lwd = 2)
       lines(c(0,no_optCopy), c(FVE[no_optCopy], FVE[no_optCopy])*100, col = "blue")
       lines(c(no_optCopy, no_optCopy), c(0, FVE[no_optCopy]*100), col = "blue")
       points(x,y,col = "green", pch = 20)
       text(no_optCopy+6, FVE[no_optCopy]*100-8, paste(" k = ", no_optCopy, ", FVE = ", round(FVE[no_optCopy]*100,digits = 3), 
            "% (threshold choice)", sep = ""))       
       no_opt = readNoPC(FVE, verbose = verbose)
       invisible(no_opt)
}