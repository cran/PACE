getScores1 = function(y, tt, mu, phi, lambda, sigma, sigma_new, noeig, error = 1, method = "CE", shrink = 0, out1, regular = 0, rho = 0){

    r1 = convertMuPhi(tt,out1, mu, phi, regular)
    muSub = r1$muSub
    phiSub = r1$phiSub
    rm(r1)

    ncohort = dim(y)[1]
    LAMBDA = diag(lambda, length(lambda), length(lambda))

    if(method == "IN"){
        xi_var = NULL
    }else{
        xi_var = vector("list",ncohort)
    }

    if(error == 1){
        
      if(sigma_new < rho){
         sigma1 = rho          
      }else{
         sigma1 = sigma_new
      }

      if(regular == 2 && method == "CE"){

          error0 = diag(sigma1, length(tt[1,]), length(tt[1,]))
          A = LAMBDA%*%t(phiSub)%*%ginv(phiSub%*%LAMBDA%*%t(phiSub)+error0)
          MU = matrix(rep(muSub, ncohort), ncohort,byrow = TRUE)
          B = y-MU
          xi_est = t(A%*%t(B))
          y_predOrig = MU+xi_est%*%t(phiSub)
          C = LAMBDA-A%*%t(LAMBDA%*%t(phiSub))
          for(i in 1:ncohort){
             xi_var[[i]] = C
          }
      }else{

          ni = sapply(1:ncohort, function(x) sum(!is.na(y[x,])))
          y_predOrig = matrix(NA, ncohort, dim(y)[2])
          xi_est = matrix(NA, ncohort, noeig)
          zeta_est = xi_est
          phii = phiSub
          mu_i = muSub
          for(i in 1:ncohort){
              if(regular != 2){
                 phii = phiSub[[i]]
                 mu_i = muSub[[i]]
              }
              yi = y[i,]
              idx = which(!is.na(yi))
              yi = yi[idx]
              if(method == "CE"){
                  error0 = diag(sigma1, length(yi), length(yi))
                  A = LAMBDA%*%t(phii)%*%ginv(phii%*%LAMBDA%*%t(phii)+error0)
                  xi_est[i,] = t(A%*%(yi-mu_i))
                  xi_var[[i]] = LAMBDA-A%*%t(LAMBDA%*%t(phii))
              }else if(method == "IN"){
                  m = length(yi)
                  ti = tt[i,]
                  ti = ti[!is.na(ti)]
                  for(k in 1:noeig){
                      prod = (yi-mu_i)*phii[,k]
                      zeta_est[i,k] = romb(ti,prod)
                      if(shrink == 0){
                         xi_est[i,k] = zeta_est[i,k]
                      }else{
                         xi_est[i,k] = lambda[k]*zeta_est[i,k]/(lambda[k]+sigma/m)
                      }
                  }
 
              }
              
              y_predOrig[i,idx] = mu_i+xi_est[i,]%*%t(phii)
          }

      }
    }else if(error == 0){
        if(regular == 2 && method == "CE"){
           
          A = LAMBDA%*%t(phiSub)%*%ginv(phiSub%*%LAMBDA%*%t(phiSub))
          MU = matrix(rep(muSub, ncohort), ncohort,byrow = TRUE)
          B = y-MU
          xi_est = t(A%*%t(B))
          y_predOrig = MU+xi_est%*%t(phiSub)
          C = LAMBDA-A%*%t(LAMBDA%*%t(phiSub))
          for(i in 1:ncohort){
             xi_var[[i]] = C
          }

        }else{

          ni = sapply(1:ncohort, function(x) sum(!is.na(y[x,])))
          y_predOrig = matrix(NA, ncohort, dim(y)[2])
          xi_est = matrix(NA, ncohort, noeig)
          phii = phiSub
          mu_i = muSub
          for(i in 1:ncohort){
              if(regular != 2){
                 phii = phiSub[[i]]
                 mu_i = muSub[[i]]
              }
              yi = y[i,]
              idx = which(!is.na(yi))
              yi = yi[idx]
              if(method == "CE"){
                  A = LAMBDA%*%t(phii)%*%ginv(phii%*%LAMBDA%*%t(phii))
                  xi_est[i,] = t(A%*%(yi-mu_i))
                  xi_var[[i]] = LAMBDA-A%*%t(LAMBDA%*%t(phii))
              }else if(method == "IN"){
                  ti = tt[i,]
                  ti = ti[!is.na(ti)]
                  for(k in 1:noeig){
                      prod = (yi-mu_i)*phii[,k]
                      xi_est[i,k] = romb(ti,prod)
                  }
 
              }
              y_predOrig[i,idx] = mu_i+xi_est[i,]%*%t(phii)
          }

        }

    }
    list(xi_est = xi_est, xi_var = xi_var, y_predOrig = y_predOrig)
}