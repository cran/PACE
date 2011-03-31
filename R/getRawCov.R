getRawCov = function(y,tt, out1new, mu, regular = 0, error = 0){

   y = as.matrix(y)
   tt = as.matrix(tt)
   out1 = unique(sort(vecData(tt)))
   ncohort = dim(y)[1]
   mu = mapX1d(out1new,mu, out1)
   indx = NULL; count = NULL;
   
   if(regular == 2){
       MU = matrix(rep(mu,each = ncohort),ncohort)    
       t1 = tt[1,]
       rm(tt)
       y = y-MU;
       cyy = t(y)%*%y/ncohort;
       cyy = as.vector(cyy)
       cxxn = cyy;
       tpairn = rbind(rep(t1, each = length(t1)), rep(t1,length(t1)));
       rm(t1,y);
 
      if(error==1){
          #tneq = apply(tpairn,2,function(x) x[1] != x[2]);
          tneq = tpairn[1,] != tpairn[2,];
          cxxn=cyy[tneq];
          tpairn=tpairn[,tneq];
      }
      win = rep(1,len = length(cxxn));
       
   }else if(regular == 1){

      MU = matrix(rep(mu,each = ncohort),ncohort)          
      ID = matrix(1, ncohort, length(out1))
      for(i in 1:ncohort){
          id = which(is.na(y[i,]))
          y[i,id] = 0
          ID[i,id] = 0
          MU[i,id] = 0
      }
      
      y = y-MU
      count = t(ID)%*%ID
      cyy = t(y)%*%y
      cyy = as.vector(cyy)
      count = as.vector(count)     
      cyy = cyy[count != 0]
      cxxn = cyy
      tpairn = rbind(rep(out1,each = length(out1)), rep(out1, length(out1)));
      rm(y, MU, out1, tt)
      tpairn = tpairn[,count != 0]
      count = count[count != 0]
      if(error == 1){
          #tneq = apply(tpairn,2,function(x) x[1] != x[2]);
          tneq = tpairn[1,] != tpairn[2,]
          cxxn=cyy[tneq];
          tpairn=tpairn[,tneq];
          count = count[tneq];
      }    
      win = rep(1,len = length(cxxn))


   }else{     
      
         xx1 = myrep(tt)
         xx2 = myrep(tt, each = FALSE)
         yy1 = myrep(y)
         yy2 = myrep(y, each = FALSE)
         rm(y)
         indx = apply(tt,1, function(x) sum(!is.na(x)))
         rm(tt)
         indx = rep(1:ncohort, indx^2)
         #id1 = myunique(xx1)$id         
         #id2 = myunique(xx2)$id
         id1 = sapply(xx1, function(x) which(out1 == x))
         id2 = sapply(xx2, function(x) which(out1 == x))
         cyy = (yy1-mu[id1])*(yy2-mu[id2])
         tpairn = rbind(xx1, xx2)
        
         if(error == 1){

             #tneq = apply(tpairn,2,function(x) x[1] != x[2]);
             tneq = tpairn[1,] != tpairn[2,]
             indx = indx[tneq]
             cxxn=cyy[tneq];
             tpairn=tpairn[,tneq];

         }else if(error == 0){
             
             cxxn = cyy

         }
         win = rep(1, len = length(cxxn))
   }
   list(tpairn = tpairn, cxxn = cxxn, indx = indx, win = win, cyy = cyy, count = count)

}
