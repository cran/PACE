gcv_mullwlsn = function(tt, ngrid = 30, regular = 0, error = 1, kernel = c("gauss","epan","rect", "quar"),rcov, verbose = "on"){

   kernel = kernel[1]
   ttt = vecData(tt)
   out1 = unique(sort(ttt))
   rm(ttt)
   a0 = min(out1)
   b0 = max(out1)
   
   h0 = getMinb(tt,out1,regular)
   if(kernel == "gauss"){
      if(is.na(h0)){
         #h0 = max(b0)  
         h0 = b0
      }
      h0 = h0*0.2
   }

   if(is.na(h0)){
       cat("Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n")
       return(list(bw_xcov = NA, gcv = NA))
   }

   rcovcount = rcov$count
   if(error == 1){
      tpairn = rcov$tpairn;
      #tneq = apply(tpairn,2,function(x) x[1] != x[2])
      tneq = tpairn[1,] != tpairn[2,]
      cyy = rcov$cyy;
      tpairn = tpairn[,tneq];
      cxxn=cyy[tneq];
      win= rep(1,len = length(cxxn))
      if(regular == 1){
          rcovcount = rcovcount[tneq]
      }
   }else{
       tpairn = rcov$tpairn
       cxxn = rcov$cxxn
       win = rcov$win
   }   
   rm(rcov,cyy,tneq)
   N = length(cxxn)
   r = b0-a0
   rm(out1)

   qq = (r/(4*h0))^(1/9);
   bw = sort(qq^(0:9)*h0);                                                                                                                    
   bw = matrix(rep(bw,2),,2)                                                                                                                        
   k0 = mykernel(0, kernel);                                                                                                                           
   out21 = seq(a0,b0,len = ngrid); 

   leave = 0
   nleave = 0
   tooSparse = 0

   while(leave == 0){
 
       gcv = rep(Inf,len = dim(bw)[1])
       for(k in 1:dim(bw)[1]){
           #cat("k = ", k, "\n")
           if(regular == 1){
                xcov = mullwlsk(bw = bw[k,], kernel = kernel, xin = tpairn, yin =cxxn , win =win , out1 = out21, out2 = out21, count = rcovcount) 
          }else{
                xcov = mullwlsk(bw = bw[k,], kernel = kernel, xin = tpairn, yin =cxxn , win =win , out1 = out21, out2 = out21) 
          }
 
          invalid = xcov$invalid
          xcov = xcov$mu
         
          if(invalid != 1){
               o21 = expand.grid(out21,out21)
               xcov = as.vector(xcov)
               #require package akima for 2-D interpolation
               #it seems to me that only linear interpolation works
               #do not allow extrapolation
               newxcov =interpp(o21[,1],o21[,2],xcov,tpairn[1,],tpairn[2,])$z
               
               rm(xcov)  
               if(regular == 1){
                   cvsum = sum((cxxn/rcovcount-newxcov)^2)
               }else{
                   cvsum = sum((cxxn-newxcov)^2)
               }
               rm(newxcov)
               bottom = 1-(1/N)*((r*k0)/bw[k])^2
             
               gcv[k] = cvsum/(bottom)^2
               tmp = gcv[gcv != Inf]
               if(length(tmp) > 1 && gcv[k] > gcv[k-1]){
                    leave = 1
                    break
               }             
          }

       }
       
       if(all(gcv == Inf)){
           if(nleave == 0 && bw[10,1] < r){
                bw_xcov = bw[10,]
                tooSparse = 1
           }else{
                cat("Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n");
                return(list(bw_xcov = NA, gcv = NA))
           }
       }else{
           bw_xcov = bw[which(gcv == min(gcv))[1],]
       }

       if(bw_xcov[1] == r){
           leave = 1
           cat("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
       }else if(bw_xcov[1] == bw[10,1] && nleave == 0){
           if((tooSparse == 1) || (sum(gcv == Inf) == 9)){
               cat("data is too sparse, retry with larger bandwidths!\n")
               h0 = bw[10,1]*1.01
           }else{
              cat("Bandwidth candidates are too small, retry with larger choices now!\n")
              h0 = bw[9,1]
           }
           newr = seq(0.5,1,by = 0.05)*r
           id = which(h0 < newr)[1]
           qq = (newr[id]/h0)^(1/9)
           bw = sort(qq^(0:9)*h0);
           bw = matrix(rep(bw,2),,2)
           if(verbose == "on"){
              cat("New bwxcov candidates:\n")
              print(bw)
           }
       }else if(bw_xcov[1] < bw[10,1] || nleave > 0){
           leave = 1
       }
       nleave = nleave+1
  }
    
    if(kernel != "gauss" && verbose == "on")
       cat(paste("GCV bandwidth choice for COV function : (", bw_xcov[1],",",bw_xcov[2],")\n", sep = ""))
    
    return(list(bw_xcov = bw_xcov, gcv = gcv))
}
