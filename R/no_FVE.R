no_FVE = function(xcov, FVE_threshold = 0.85, disp = FALSE){
   d = eigen(xcov, only.values = FALSE)$values
   idx = which(Im(d) != 0)      #find indices for imaginary eigenvalues
   if(length(idx) >  0){
      stop(paste(length(idx), "eigenvalues are complex. The estimated auto-covariance surface is not symmetric!"));
   } 
   idx = which(d <= 0)
   if(length(idx) > 0){
       if(disp){
         cat(paste(length(idx), "real eigenvalues are negative or zero and are removed!\n"))
       }
       d = d[d > 0]
   }
   lambda = sort(d, decreasing = TRUE)   
   FVE = cumsum(lambda)/sum(lambda)
   no_opt = which(FVE > FVE_threshold)[1]
   list(no_opt = no_opt, FVE = FVE, lambda = lambda)
}