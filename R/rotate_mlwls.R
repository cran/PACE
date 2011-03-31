rotate_mlwls = function(bw, kernel =c("gauss","epan","rect", "quar"), xin,yin, win, d, npoly = 1){
   active = win != 0
   xin = xin[,active]
   yin = yin[active]
   win = win[active]

   #rotating coordinates of predictors by pi/4
   R = sqrt(2)/2*matrix(c(1,-1,1,1),2,2,byrow = TRUE)
   xn = R%*%xin
   yn = yin
   dn = R%*%d
   
   invalid = 0
   #minimizing local weighted least squares
   m = dim(d)[2]
   mu = numeric(m)
   for(i in 1:m){
      #locating local window
      if(kernel != "gauss" && kernel != "gausvar"){
          list1 = which(xn[1,] >= dn[1,i]-bw[1] & xn[1,] <= dn[1,i]+bw[1])
          list2 = which(xn[2,] >= dn[2,i]-bw[2] & xn[2,] <= dn[2,i]+bw[2])
          ind = intersect(list1,list2)
      }else{
          ind = 1:dim(xn)[2]
      }

      lx = xn[,ind]
      ly = yn[ind]
      lw = win[ind]
      #computing weight matrix
      if(length(ly) >= npoly+1){
           llx = rbind((lx[1,]-dn[1,i])/bw[1],(lx[2,]-dn[2,i])/bw[2])
           #deciding the kernel to use
           if(kernel == "epan"){
              w = lw*(1-llx[1,]^2)*(1-llx[2,]^2)*(9/16)
           }else if(kernel == "rect"){
              w = lw*rep(1,dim(lx)[2])/4
           }else if(kernel == "gauss"){
              w = lw*dnorm(llx[1,])*dnorm(llx[2,])
           }else if(kernel == "gausvar"){
              w = lw*dnorm(llx[1,])*(1.25-0.25*llx[1,]^2)*dnorm(llx[2,])*(1.5-0.5*llx[2,]^2)
           }else if(kernel == "quar"){
              w = lw*((1-llx[1,]^2)^2)*((1-llx[2,]^2)^2)*(225/256)
           }
           W = diag(w, length(w), length(w))
           #computing design matrix
           X = matrix(1,length(ly),3)
           #X[,1] = rep(1,len = length(ly))
           X[,2] = (lx[1,]-dn[1,i])^2
           X[,3] = lx[2,]-dn[2,i]
    
           beta = ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%ly
           mu[i] = beta[1]


      }else if(length(ly) == 1){
           mu[i] = ly
      }else{
           cat("No points in local window, please increase bandwidth!\n")
           invalid = 1;
           return(list(invalid = invalid, mu = NULL))
      }

   }
   list(invalid = invalid, mu = mu)
}
