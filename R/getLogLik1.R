getLogLik1 = function(y,lambda,muSub, phiSub, sigma, regular = 0){
 
    invalid = 0
    logLik1 = 0
    k = length(lambda)
   
    LAMBDA = diag(lambda, length(lambda), length(lambda))
    
    if(regular == 2){
        mm = dim(phiSub)[1]
        phiSub = phiSub[,1:k]      
        if(!is.matrix(phiSub))
           phiSub = matrix(phiSub,mm,k)
        
        error0 = diag(sigma, length(muSub), length(muSub))
        Sigma_y = phiSub%*%LAMBDA%*%t(phiSub)+error0
        detSigma_y = det(Sigma_y)
        if(detSigma_y < 10^(-10)){
           invalid = 1
           logLik1 = NULL
           return(list(invalid = invalid, logLik1 = logLik1))
        }
        for(i in 1:dim(y)[1]){
          yi = y[i,]
          yi = yi[!is.na(yi)]
          logLik1 = log(detSigma_y)+(yi-muSub)%*%ginv(Sigma_y)%*%(yi-muSub)+logLik1
        }
    }else{
        for(i in 1:dim(y)[1]){
          phii = phiSub[[i]]
          mm = dim(phii)[1]
          
          phii = phii[,1:k]
          if(!is.matrix(phii))
            phii = matrix(phii, mm,k)
          
          mu_i = muSub[[i]]
          error0 = diag(sigma, length(mu_i), length(mu_i))
         
          Sigma_y = phii%*%LAMBDA%*%t(phii)+error0
          detSigma_y = det(Sigma_y)
    
          if(detSigma_y < 10^(-10)){
             
             invalid = 1
             logLik1 = NULL
             return(list(invalid = invalid, logLik1 = logLik1))
          }else{
             yi = y[i,]
             yi = yi[!is.na(yi)]
             logLik1 = logLik1 + log(detSigma_y)+(yi-mu_i)%*%ginv(Sigma_y)%*%(yi-mu_i)
          }

        }
    }
    return(list(invalid = invalid, logLik1 = as.numeric(logLik1)))
}
