lwls = function(bw, kernel = c("gauss","epan","rect", "quar"), nwe = 0, 
npoly = nder+1, nder = 0, xin, yin, win, xou, bwmuLocal = 0){

   if(npoly < nder)
      stop("Degree of polynomial should be no less than the order of derivative!")
   
   kernel = kernel[1]
   require(MASS)
   actobs = which(win != 0)
   xin = xin[actobs]
   yin = yin[actobs]
   win = win[actobs]
   invalid = 0

   aa = 1
   if(nwe == -1)
      aa = 4
   mu = numeric(length(xou))
   gap = mu

   search_ht = function(t1,h0){
     tmp = -(t1-75)^2/6000
     h0*(1+exp(tmp))

   }

   if(bw > 0){
      if(bwmuLocal){
         bw = search_ht(xou,bw)
      }else{
         bw = rep(bw,length(xou)) 
      }
   }else{
      stop("Bandwidth choice for mu(t) and/or its derivative must be positive!")
   }
   
   #LWLS with different weight functions
   for(i in 1:length(xou)){

       #(3-1) Locating local window
       if(kernel != "gauss" && kernel != "gausvar"){
          idx = xin <= xou[i] + aa*bw[i] & xin >= xou[i]-aa*bw[i]       
       }else{
          idx = 1:length(xin)
       }
       lx = xin[idx]
       ly = yin[idx]
       lw = win[idx]

       if(length(unique(lx)) >= (npoly+1)){

          #Sepcify weight matrix
          llx = (lx-xou[i])/bw[i]
          
          if(kernel == "epan"){
              w = lw*(1-llx^2)*0.75
          }else if(kernel == "rect"){
              w = lw
          }else if(kernel == "optp"){
              w = lw*(1-llx^2)^(nwe-1)
          }else if(kernel == "gauss"){
              w = lw*dnorm(llx)
          }else if(kernel == "gausvar"){
              w = lw*dnorm(llx)*(1.25-0.25*llx^2)
          }else if(kernel == "quar"){
              w = lw*(15/16)*(1-llx^2)^2
          }else{
              cat("Invalid kernel, Epanechnikov kernel is used!\n")
              w = lw*(1-llx^2)*0.75
          }
          W = diag(w, length(w), length(w))
          # Define design matrix
          dx = matrix(1,length(lx),npoly+1)
          for(j in 1:npoly){
            dx[,j+1] = (xou[i]-lx)^j
          }       
          
          p = ginv(t(dx)%*%W%*%dx)%*%t(dx)%*%W%*%ly  

          #Find estimate
          mu[i] = p[(nder+1)*gamma(nder+1)*((-1)^nder)]          
          
       }else{
          gap[i] = 1
          invalid = 1
       }

   }
   indx = which(gap == 0)
   if((length(indx)>= 0.9*length(xou))&& (length(indx)<length(xou))){
        mu1 = mu[indx]
        rr = myunique(xou[indx])
        mu=interp1(rr$out1,mu1[rr$id],xou);
   }else if(length(indx) < 0.9*length(xou)){
        mu = NULL
        cat("Too many gaps, please increase bandwidth!\n")
        invalid = 1
   }
   return(list(invalid = invalid, mu = mu))
}
