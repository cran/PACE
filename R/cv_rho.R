cv_rho = function(y, tt, mu, phi, lambda, sigma, sig1, noeig, error = 1, method = "CE", shrink = 0, out1, regular = 0, rho = "cv", ni, tjID, verbose = "on"){

   subID = which(ni > 1)  #find indices for subjects with ni >= 2
  
   r1 = convertMuPhi(tt,out1, mu, phi, regular)
   muSub = r1$muSub
   phiSub = r1$phiSub
   rm(r1)

   LAMBDA = diag(lambda, length(lambda), length(lambda))
   cv = rep(Inf, length(rho))

   for(k in 1:length(rho)){
      cv[k] = getScores2(y, tt, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, muSub, phiSub, LAMBDA, rho[k], subID, tjID)     
   }
  
   IDopt = which(cv == min(cv))[1]
   rho_opt = rho[IDopt]
   if(verbose == "on"){
      cat("Best rho from CV method:", rho_opt, "with cv =", cv[IDopt],"\n")
   }
   rho_opt

}
