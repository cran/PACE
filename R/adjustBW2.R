adjustBW2 = function(kernel=c("gauss","epan","rect", "quar"), 
bw_xcov, npoly = nder+1 , nder = 0, regular = 0, verbose = "on")
{
    kernel = kernel[1]
    if(kernel == "gauss"){
        if(regular == 2){
          bwxcov_fac = c(1.1,0.8,0.8)
        }else{
          bwxcov_fac = c(1.1, 1.2, 2)
        }
       if(nder > 2){
          facID = 3
       }else if(nder >= 0 && nder <= 2){
          facID = nder+1
       }else{
          facID = 1
       }
       bw_xcov = bw_xcov*bwxcov_fac[facID]
       if(verbose == "on")
         cat("Adjusted GCV bandwidth choice for COV function (npoly = ", npoly, "): (", bw_xcov[1],",", bw_xcov[2], ")\n")
    }
    bw_xcov   

}
