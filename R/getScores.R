getScores = function(y, tt, mu, phi, lambda, sigma, noeig, error =  1, 
                     method = "CE", shrink = 0, out1, regular = 0, rho, verbose = "on")
{
  
   n = dim(y)[1]
   sig1 = sigma

   #Compute iterative residuals
   if(rho != -1){
      for(j in 1:2){
         SE_i = rep(Inf, len = n)         
         yo = getOriCurves(y,tt,mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular)$y_predOrig
         for(i in 1:n){
            yi = y[i,]
            idx = which(!is.na(yi))
            yi = yi[idx]
            yo_i = yo[i,idx]
            SE_i[i] = mean((yi-yo_i)^2)
         }
         sig1 = mean(SE_i)
      }
   }
   
   if(is.numeric(rho)){
       if(rho >= 0 || rho == -1){
           rho_opt = rho
       }else{
           cat("rho should not be negative! Reset it to cv choice now!\n")
           rho = "cv"
       }
   } 

#   if(rho == "cv"){
    if(is.character(rho) && length(grep("cv", rho)) > 0){
       #Compute gamma
       TT = range(tt, na.rm = TRUE)
       TT = TT[2]-TT[1]
       gamma = ((romb(out1, mu^2)+ sum(lambda))/TT)^(0.5)
       alpha = seq(0.01, 0.22, len = 50)
       rho = gamma*alpha

       #ni = numeric(n)
       #tjID = numeric(n)
       #for(i in 1:n){
       #   ni[i] = sum(!is.na(tt[i,]))
       #   if(ni[i] > 1)
       #      tjID[i] = sample(1:ni[i],1)
       #}

      #default is "cv", non-randomized leave-one-measurement-out
      isRandom = 0;
      if(rho == "cv-random"){
         isRandom = 1;   #randomized leave-one-measurement-out
      }
     
      rr = getTimeID(tt, n, isRandom)
      tjID = rr$tjID
      ni = rr$ni
      rm(rr)     
 
      #find optimal rho from subjects with ni >= 2
      if(all(ni == 1) || error == 0){
          rho_opt = min(rho)
      }else{
         
          rho_opt = cv_rho(y,tt, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, rho, ni, tjID, verbose)
      }
   }
   
   r = getScores1(y,tt,mu,phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, rho_opt)
   list(xi_est = r$xi_est, xi_var = r$xi_var, y_predOrig = r$y_predOrig, rho_opt = rho_opt, sig1 = sig1)
}
