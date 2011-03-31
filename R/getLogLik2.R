getLogLik2 = function(y,tt, sigma, noeig, method = "CE", shrink = 0, regular = 0, muSub, phiSub, lambda){

  k = length(lambda)
  LAMBDA = diag(lambda, length(lambda), length(lambda))
  
  if(regular == 2 && method == "CE"){

      mm = dim(phiSub)[1]
      phiSub = phiSub[,1:k]
      if(!is.matrix(phiSub))
         phiSub = matrix(phiSub,mm,k)
      error0 = diag(sigma, length(tt[1,]), length(tt[1,]))
      A = LAMBDA%*%t(phiSub)%*%ginv(phiSub%*%LAMBDA%*%t(phiSub)+error0)
      MU = matrix(rep(muSub, dim(y)[1]), dim(y)[1],byrow = TRUE)
      B = y-MU
      xi_est = t(A%*%t(B))
      y_pred = MU+xi_est%*%t(phiSub)
      res = y-y_pred
      res = (res^2)/sigma
      logLik2 = sum(res)
  }else{
      xi_est = numeric(noeig)
      zeta_est = xi_est
      logLik2 = 0
      phii = phiSub
      mu_i = muSub
      for(i in 1:dim(y)[1]){
         if(regular != 2){
            phii = phiSub[[i]]
            mu_i = muSub[[i]]
         }

         mm = dim(phii)[1]
         phii = phii[,1:k]
         if(!is.matrix(phii))
            phii = matrix(phii, mm,k)
     
         yi = y[i,]
         yi = yi[!is.na(yi)]
         
         if(method == "CE"){
            error0 = diag(sigma, length(yi), length(yi))
            A = LAMBDA%*%t(phii)%*%ginv(phii%*%LAMBDA%*%t(phii)+error0)
            xi_est = t(A%*%(yi-mu_i))
         }else if(method == "IN"){
            m = length(yi)
            ti = tt[i,]
            ti = ti[!is.na(ti)]
            for(k in 1:noeig){
               prod = (yi-mu_i)*phii[,k]
               zeta_est[k] = romb(ti,prod)
               xi_est[k] = lambda[k]*zeta_est[k]/(lambda[k]+sigma/m)
               if(shrink == 0)
                  xi_est[k] = zeta_est[k]
            }
         }
         y_pred = mu_i+xi_est%*%t(phii)
         logLik2 = sum((yi-y_pred)^2)/sigma+logLik2
      }
  }

  as.numeric(logLik2)
}
