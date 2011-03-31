gcv_lwls = function(yy, tt, kernel = c("gauss","epan","rect", "quar"), nwe = 0,
npoly = nder+1, nder = 0, regular = 0, verbose = "on", bwmuLocal = 0){
    kernel = kernel[1]
    r = range(tt, na.rm = TRUE)
    r = r[2]-r[1]
    N = length(tt)
    if(regular == 0){
       dstar = minb(tt,2+npoly)
       if(dstar > r/4){
          dstar = dstar*.75
          cat("The min bandwidth choice is too big, reduce to ", dstar, " now!\n");
       }
       h0 = 2.5*dstar
    }else if (regular == 1){
       h0 = minb(tt,1+npoly)*2
    }else{
       h0 = minb(tt,1+npoly)*1.5
    }

    if(is.na(h0)){
        if(kernel == "gauss"){
           h0 = 0.2*r;
        }else{
           cat("Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
           return(list(bopt = NA, gcv = NA))
        }
    }
 
    h0 = min(h0, r)
    qq = (r/(4*h0))^(1/9)
    #out1 = unique(tt)
    #idx = sapply(tt, function(x) which(out1 == x))    
    tmp = myunique(tt)
    out1 = tmp$out1
    idx = tmp$id    
    rm(tmp)

    exps = 0:9
    bw = sort(qq^exps*h0)
   
    win1 = rep(1,N)
    k0 = mykernel(0,kernel)
    
    leave = 0
    nleave = 0
    tooSparse = 0
    while(leave == 0){
        gcv = rep(Inf,len = length(bw))
        for(k in 1:length(bw)){
            #print(paste("k = ", k))
            if(length(out1) > 101){
                out21 = seq(min(out1),max(out1),len = 101)
                res = lwls(bw[k], kernel, nwe, npoly, nder, tt, yy,win1,out21,bwmuLocal)
                if(res$invalid == 0)
                    newmu = interp1(out21,res$mu,tt)
            }else{
                res = lwls(bw[k], kernel, nwe, npoly, nder, tt, yy,win1,out1,bwmuLocal)
                if(res$invalid == 0)
                    newmu = res$mu[idx]
            }
            if(res$invalid == 0){
                cvsum = as.numeric((yy-newmu)%*%(yy-newmu))
                gcv[k] = cvsum/(1-(r*k0)/(N*bw[k]))^2
                if(k > 1 && (gcv[k] > gcv[k-1])){
                   leave = 1
                   break
                }
            }
        }
       
        if(all(gcv == Inf)){
            if(nleave == 0 && bw[10] < r){
                bopt = bw[10]
                tooSparse = 1;
            }else{
                cat("Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
                return(list(bopt = NA, gcv = NA))               
            }            
        }else{
            bopt = bw[which(gcv == min(gcv))[1]]
        }

        if(bopt == r){
           leave = 1
           cat("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
        }else if(bopt == bw[10] && nleave == 0){

           if((tooSparse == 1)||(sum(gcv == Inf) == 9)){
               cat("data is too sparse, retry with larger bandwidths!\n")
               h0 = bw[10]*1.01
           }else{
              cat("Bandwidth candidates are too small, retry with larger choices now!\n")
              h0 = bw[9]
           }
           newr = seq(0.5,1,by = 0.05)*r
           id = which(h0 < newr)[1]
           qq = (newr[id]/h0)^(1/9)
           exps = 0:9
           bw = qq^exps*h0
           bw = sort(bw)
           if(verbose == "on"){
              cat("New bwmu candidates:\n")
              print(bw)
           }
        }else if(bopt < bw[10] || nleave > 0){
           leave = 1
        }
        nleave = nleave +1 
    }
    if(kernel != "gauss" && verbose == "on"){
        cat("GCV choice for mean function (npoly = ", npoly, "): ", bopt, "\n")
    }
    list(bopt = bopt, gcv = gcv)
}
