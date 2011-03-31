getBINnum = function(n,m,regular = 0, verbose = "on"){

   numBin = NULL;
   if(m<=20){
      if(regular == 0){
         str = "Median of ni";
      }else{
         str = "Maximum of ni";
      }
      if(verbose== "on")
          cat(str, "is no more than 20!No binning is performed!\n");
      return(numBin);
   }
   if(m > 400)
      numBin = 400;
    
   if(n > 5000){
      mstar = max(20,(5000-n)*19/2250+400);
      if(m > mstar)
         numBin = ceiling(mstar)
   }

   if(is.null(numBin) && verbose == "on")
      cat("No binning is needed.\n")
   return(numBin);
}
