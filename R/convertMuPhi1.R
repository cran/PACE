convertMuPhi1 = function(tt, out1, mu, phi, regular = 0){

      if(regular == 2){
         tmp = tt[1,];
         tmp = tmp[!is.na(tmp)]
         muSub = interp1(out1,mu,tmp)
         phiSub = interp11(out1,phi,tmp)
         if(!is.matrix(phiSub)){
            phiSub = matrix(phiSub, 1,length(phiSub));
         }
      }else{
         muSub = vector("list",dim(tt)[1])
         phiSub = muSub
         for(i in 1:dim(tt)[1]){
           tmp = tt[i,]
           tmp = tmp[!is.na(tmp)]
           muSub[[i]] = interp1(out1,mu,tmp)
           phiSub[[i]] = interp11(out1,phi,tmp)
           if(!is.matrix(phiSub[[i]])){
              phiSub[[i]] = matrix(phiSub[[i]], 1,length(phiSub[[i]]));
           }
         }
      }
      list(muSub = muSub, phiSub = phiSub)
}
   
