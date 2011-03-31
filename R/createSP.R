createSP = function(FVE, no_opt, yname = "y"){

   x11()
   x = 0:length(FVE)
   y = c(0,FVE)*100
   plot(x,y, type = "p", lwd = 3, xlim = c(0,length(FVE)+1), ylim = c(0,101),
        main = paste("Fraction of variance explained by No. of PC for function ", yname, sep = ""),
        xlab = "No. of Principal Components", ylab = "FVE(%)", font.lab = 2)
   lines(x,y,lty = 2, col = "red", lwd = 2)
   lines(c(0,no_opt),c(FVE[no_opt], FVE[no_opt])*100,col = "blue")
   lines(c(no_opt, no_opt),c(0,FVE[no_opt])*100 ,col = "blue")
   points(x,y,col = "green", pch = 20)
   text(no_opt+6, FVE[no_opt]*100-8, paste(" k = ", no_opt, ", FVE = ", round(FVE[no_opt]*100, digits = 3), 
        "% (final choice)", sep = ""))

}