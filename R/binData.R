binData = function(y, tt, regular = 0, verbose = "on", numBins = NULL){

   n = nrow(tt);
   ni = apply(tt,1,function(x) sum(!is.na(x)))

   if(regular == 0){
      m = median(ni);
   }else{
      m = max(ni);
   }
   if(is.null(numBins)){
      numBins = getBINnum(n,m,regular, verbose);
   }else if(numBins <= 0){
      cat("Number of bins must be positive integer! Reset to default number of bins now!\n")
      numBins = getBINnum(n,m,regular, verbose);
   }

   if(is.null(numBins)){
      newy = NULL;
      newt = NULL;
      return(list(newy = newy, newt = newt))
   }

   numBins = ceiling(numBins);
   newt = matrix(NA, n, numBins)
   newy = matrix(NA, n, numBins)
   
   a0 = min(tt, na.rm = TRUE)
   b0 = max(tt, na.rm = TRUE)
  
   if(verbose == "on")
      cat("Start binning: \n")
   for(i in 1:n){
       res = binning(tt[i,], y[i,], numBins, 1,1,a0,b0)
       npts = sum(!is.na(res$midpoint))
       newt[i,1:npts] = res[["midpoint"]][1:npts]
       newy[i,1:npts] = res[["newy"]][1:npts]
   }
   if(verbose == "on")
      cat("Number of bins selected for the input data in the time domain [", a0, ",",b0, "] is ", numBins, " .\n")
   return(list(newy = newy, newt = newt))
}