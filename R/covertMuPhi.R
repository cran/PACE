convertMuPhi = function(tt, out1, mu, phi, regular = 0){
  
   if(regular == 2){
      muSub = mapX1d(out1,mu,tt[1,])
      phiSub = mapX1d(out1,phi,tt[1,])
      if(!is.matrix(phiSub)){
              phiSub = matrix(phiSub, 1,length(phiSub));
      }
   }else{
      muSub = vector("list",dim(tt)[1])
      phiSub = muSub
      for(i in 1:dim(tt)[1]){
        tmp = tt[i,]
        tmp = tmp[!is.na(tmp)]
        muSub[[i]] = mapX1d(out1,mu,tmp)
        phiSub[[i]] = mapX1d(out1,phi,tmp)
        if(!is.matrix(phiSub[[i]])){
              phiSub[[i]] = matrix(phiSub[[i]], 1,length(phiSub[[i]]));
        }

      }
   }
   list(muSub = muSub, phiSub = phiSub)
}
