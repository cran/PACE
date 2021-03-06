no_AIC = function(y,tt, mu, bw_xcov, ngrid = 51, regular = 0, maxk = 20, AB_method = "full", 
                 method = "CE", shrink = 0, out1, out21, kernel = c("gauss","epan","rect", "quar"),
                 error = 1, cut = 1, rcov, xcov, npoly =1){
     aic = rep(Inf, maxk)
     if(error == 1){
         r = pc_covE(tt, out1, bw_xcov, ngrid, cut, kernel, rcov, npoly)
         invalid = r$invalid
         sigma = r$sigma
     }else{
         sigma = 0
     }
     
    r = getEigens(xcov, out1, out21, maxk)
    lambda = r$lambda
    phi = r$phi
    eigens = r$eigen
    noeig = r$noeig
    rm(r)

    r1 = convertMuPhi(tt,out1, mu, phi, regular)
    muSub = r1$muSub
    phiSub = r1$phiSub
    rm(r1)

    if(noeig < maxk){
        cat(paste("Max number of PC for AIC selection is no more than", noeig, "! You may want to increase maxk = ", maxk, "to a larger number to include greater flexibility.\n"))
        maxk = noeig 
    }
    
    for(k in 1:maxk){
       sig = lambda[1:k]
       if(AB_method == "full"){
            r = getLogLik1(y,sig, muSub, phiSub, sigma, regular)
            invalid = r$invalid
            logLik = r$logLik
            rm(r)
            if(invalid){
                cat("The covariance matrix of the estimated function is nearly singular! Reset to FVE method now!\n")
                return(list(no_opt = NULL, aic = NULL))
            }
       }else{
            logLik = getLogLik2(y,tt,sigma,k,method,shrink,regular,muSub, phiSub, sig)
       }
       aic[k] = logLik + 2*k
       tmp = aic[aic != Inf]
       if(length(tmp) > 1 && k > 1){
          if(aic[k] > aic[k-1])
             break
       }
    }
    aic = aic[aic != Inf]
    no_opt = which(aic == min(aic))[1]
    list(no_opt = no_opt, aic = aic)   
} 