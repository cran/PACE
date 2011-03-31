cv_mullwlsn = function(y,tt,mu,ntest = 30,ngrid = 30,regular = 0,error = 1,
                       kernel = c("gauss","epan","rect", "quar"),rcov, verbose = "on"){
  
   kernel = kernel[1]
   ncohort = dim(tt)[1]
   ind = apply(tt,1,function(x) sum(!is.na(x)))
   ind = rep(1:length(ind), times = ind)
   ttt = vecData(tt)
   #yy = vecData(y)
   lint = range(ttt)
   lint = lint[2]-lint[1]
   win1 = rep(1, len = length(ttt))
   out1 = unique(sort(ttt))
   n = length(out1)
   rm(ttt)

   test = sort(sample(1:ncohort, ntest))
   ni = sapply(test, function(x) sum(!is.na(y[x,])))
   act = sum(ni > 1)
   ni = sapply(1:ncohort, function(x) sum(!is.na(y[x,])))

   if(error == 1){
      tpairn = rcov$tpairn;
      #tneq = apply(tpairn,2,function(x) x[1] != x[2])
      tneq = tpairn[1,] != tpairn[2,]
      cyy = rcov$cyy;
      rcov$tpairn = tpairn[,tneq];
      rcov$cxxn=cyy[tneq];
      rcov$win= rep(1,len = length(rcov$cxxn))
      if(regular == 1){
          rcov$count = rcov$count[tneq]
      }
      rm(tpairn, tneq, cyy)
   }
   
   dstar = min(lint/4, getMinb(tt,out1,regular))

   #create 10 bandwidth candidates 
   nbw = 11
   bw = numeric(nbw-1)
   for(i in 1:(nbw-1)){
     bw[i] = 2.5*lint/n*(n/2.5)^((i-1)/(nbw-1))
   }

   bw = bw-min(bw)+dstar
   if(kernel == "gauss")
     bw = bw*0.5
   bw = matrix(rep(bw,2),,2)
   
   out21 = seq(min(out1), max(out1), by = lint/ngrid)
   out22 = out21
 
   crossv = numeric(nbw-1)
   count = numeric(nbw-1)
   
   enter = max(diff(out1))
   begin = which(bw[,1] >= enter)
   i = begin[1]
   continu = 1
   #browser()
   while(continu == 1){
        for(j in 1:ntest){
            omit = test[j]
            tomit = tt[omit,]
            yomit = y[omit, !is.na(tomit)]
            tomit = tomit[!is.na(tomit)]

       
            if ((error == 1 && ni[omit] > 1) || error == 0){
              
                if(regular == 0){
                     winomit = rcov$win
                     winomit[rcov$indx == omit] = 0
                     xcovomit = mullwlsk(bw = bw[i,],kernel = kernel,xin = rcov$tpairn,yin =rcov$cxxn , 
                                         win =winomit , out1 = tomit, out2 = tomit) 
                 }else if(regular == 2){
                     rcovnew = getRawCov(y[-omit,], tt[-omit,], out1, mu, regular, error)
                     winomit = rcovnew$win
                     xcovomit =  mullwlsk(bw = bw[i,],kernel = kernel,xin = rcovnew$tpairn,yin =rcovnew$cxxn ,
                                         win =winomit , out1 = tomit, out2 = tomit)                     
                }else if(regular == 1){
                     rcovnew = getRawCov(y[-omit,], tt[-omit,], out1, mu, regular, error)
                     winomit = rcovnew$win
                     xcovomit =  mullwlsk(bw = bw[i,],kernel = kernel,xin = rcovnew$tpairn,yin =rcovnew$cxxn ,
                                         win =winomit , out1 = tomit, out2 = tomit, rcovnew$count)                     
                }   
                invalid = xcovomit$invalid
                xcovomit = xcovomit$mu 

                if(invalid == 0){
                 
                   pred = numeric(0)
                   obs = numeric(0)
                   k = 1
                   for(m1 in 1:length(tomit)){
                       for(m2 in 1:length(tomit)){
                            if((error == 1 && m1 != m2) || error == 0){
                                pred = c(pred, xcovomit[m1,m2])
                                obs[k] =  (yomit[m1]-mu[out1 == tomit[m1]])*(yomit[m2]-mu[out1 == tomit[m2]])
                 
                                k = k+1         
                            }
                       }
                   }
 
                   crossv[i] = crossv[i]+sum((obs-pred)^2)
                   count[i] = count[i]+1
                   
                }
           }

        }
        print(crossv)
        if(i == begin[1])
            continu == 1
        else if(error == 1 && count[i-1]/act < 0.95)
            continu = 1
        else if(error == 0 && count[i-1]/ntest < 0.95)
            continu = 1
        else
            continu = (i < nbw-1) && (crossv[i] < crossv[i-1])
        if(continu == 1)
           i = i+1 
       
   }
   opt = i-1
   bw_xcov = bw[opt,]
   if(verbose == "on")
       cat(paste("CV bandwidth choice for COV function: (", bw_xcov[1], ",", bw_xcov[2], ")\n", sep = "")) 
   bw_xcov
}
