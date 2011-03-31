cvfda_lwls = function(x,tt,kernel = c("gauss","epan","rect", "quar"), nwe = 0, 
                      npoly = nder+1, nder = 0, regular = 0, verbose = "on", bwmuLocal = 0)
{
  x = as.matrix(x)
  tt = as.matrix(tt)
  kernel = kernel[1]
  ncohort = nrow(tt)
  
  ttt = vecData(tt)      #vectorized time points
  xx = vecData(x)        #vectorized repeated measurements

  ind = apply(tt,1,function(x) sum(!is.na(x)))  
  ind = rep(1:length(ind), times = ind)   

  tttn = ttt
  xxn = xx

  rt = range(ttt)
  rang = rt[2]-rt[1]
  dstar = minb(ttt,npoly+2);

  if(regular != 2){
     h0 = 2.5*dstar
  }else{
     h0 = dstar
  }
  
  if(h0 > rang/4){
     h0 = h0 * .75
     cat("The min bandwidth choice is too big, reduce to ",h0, "!\n")
  }
  
  #create 10 bandwidth candidates
  nbw = 11
  bw = numeric(nbw-1)
  n = length(unique(ttt))
  for(i in 1:(nbw-1)){
     bw[i] = 2.5*rang/n*(n/2.5)^((i-1)/(nbw-1))
  }
  bw = bw-min(bw)+h0
  
  if(regular == 2){
      aves = numeric(dim(tt)[2])
      for(i in 1:ncohort){
         aves = aves + tt[i,]/ncohort
      }
  }

  cv = numeric(length(bw)-1)
  count = cv

  for(j in 1:(nbw-1)){
     #cat("j = ", j,"\n")
     for(i in 1:ncohort){
         #cat("i = ", i,"\n")
         out = ttt[ind == i]
         obs = xx[ind == i]
         win = rep(1,len = length(ttt))
         win[ind == i] = 0    
         if(regular == 2){
            xxn = (aves*ncohort-tt[i,])/(ncohort-1)
            tttn = tt[1,]
            win  = rep(1,len = length(tt[1,]))      
         }

         
         res = lwls(bw[j],kernel,nwe,npoly,nder,tttn,xxn,win, out, bwmuLocal)

         if(res$invalid == 0){
             cv[j] = cv[j] + sum((obs-res$mu)^2)
             count[j] = count[j] +1
         }
     }
     
  }
  cv = cv[which(count/ncohort > 0.9)]
  bw = bw[which(count/ncohort > 0.9)]
  bopt = bw[which(cv == min(cv))[1]]
  if(verbose == "on")
      cat("CV bandwidth choice of mean function: ", bopt, "\n")
  bopt

}
