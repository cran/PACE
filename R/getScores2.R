getScores2 = function(y, tt, mu, phi, lambda, sigma, sig1, noeig, error = 1, method = "CE", shrink = 0, out1, 
                      regular = 0, muSub, phiSub, LAMBDA, rho = "cv", subID, tjID)
{

          cv = 0;
          ncohort = length(subID)  #here ncohort is the number of subjects with ni >= 2 only
          tjID = tjID[subID]

          if(sig1 < rho){
             sigma1 = rho
          }else{
             sigma1 = sig1
          }   
         
          #ni = sapply(subID, function(x) sum(!is.na(y[x,])))

          xi_est = matrix(NA, ncohort, noeig)
          zeta_est = xi_est
          phii = phiSub
          mu_i = muSub
          for(i in 1:ncohort){
              if(regular != 2){
                 phii = phiSub[[subID[i]]]
                 mu_i = muSub[[subID[i]]]
              }else{
                 phii = phiSub
                 mu_i = muSub
              }
                      
              yi = y[subID[i],]
              idx = which(!is.na(yi))
              yi = yi[idx]
              ti = tt[subID[i],idx]
              yij = yi[tjID[i]]
              tij = ti[tjID[i]]
              
              yi = yi[-tjID[i]]
              ti = ti[-tjID[i]]
              mu_ij = mu_i[tjID[i]]
              mu_i =  mu_i[-tjID[i]]
              
              if(is.matrix(phii)){
                  phii_j = phii[tjID[i],]
                  phii = phii[-tjID[i],]
              }else{
                  phii_j = phii[tjID[i]]
                  phii = phii[-tjID[i]]
              }

              if(!is.matrix(phii)){
                 phii = matrix(phii, length(yi), dim(LAMBDA)[1])
              }
              
              if(method == "CE"){
                  error0 = diag(sigma1, length(yi), length(yi))
                  A = LAMBDA%*%t(phii)%*%ginv(phii%*%LAMBDA%*%t(phii)+error0)
                  xi_est[i,] = t(A%*%(yi-mu_i))
              }else if(method == "IN"){
                  m = length(yi)
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

              cv = cv+(mu_ij+xi_est[i,]%*%phii_j-yij)^2
             
          }
          cv
}