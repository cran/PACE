mullwlsk = function(bw,kernel = c("gauss","epan","rect", "quar"),xin,yin,win = rep(1, length(xin)),out1,out2,count = NULL){
    require(MASS) #for generalized inverse
    kernel = kernel[1]
    active = which(win != 0)
    xin = xin[,active]
    yin = yin[active]
    win = win[active]
    invalid = 0
    mu = matrix(NA, length(out2), length(out1))
    for(i in 1:length(out2)){
#        cat("i = ", i, "\n")
        for(j in i:length(out1)){
#           cat("j = ", j ,"\n")
           #locating local window
           if(kernel != "gauss"){
               list1 = (xin[1,] >= out1[j]-bw[1]-10^(-6)) & (xin[1,] <= out1[j] + bw[1]+10^(-6))
               list2 = (xin[2,] >= out2[i]-bw[2]-10^(-6)) & (xin[2,] <= out2[i] + bw[2]+10^(-6))
               #ind = which(list1 == TRUE & list2 == TRUE)
               ind = list1 & list2 
           }else{
               #ind = 1:dim(xin)[2]
               ind = !logical(dim(xin)[2])
           }
           lx = xin[,ind]
           ly = yin[ind]
           lw = win[ind]
           #computing weight matrix
           if(dim(unique(t(lx)))[1]>=3){     #at least 3 unique number of time pairs in the local window
               llx = rbind((lx[1,]-out1[j])/bw[1], (lx[2,]-out2[i])/bw[2])
               #deciding the kernel used
               k = dim(llx)[2]
               #indd = 1:k
               if(kernel == "epan"){
                   temp = lw*(1-llx[1,]^2)*(1-llx[2,]^2)*(9/16)               
               }else if(kernel == "rect"){
                   temp = lw*rep(1,dim(lx)[2])/4
               }else if(kernel == "gauss"){
                   temp = lw*dnorm(llx[1,])*dnorm(llx[2,])
               }else if(kernel == "gausvar"){
                   temp = lw*dnorm(llx[1,])*(1.25-0.25*llx[1,]^2)*(dnorm(llx[2,])*(1.5-0.5*llx[2,]^2))
               }else if(kernel == "quar"){
                   temp = lw*((1-llx[1,]^2)^2)*((1-llx[2,]^2)^2)*(225/256)
               }
 
               W = diag(temp, length(temp),length(temp))
  
               #computing design matrix
               X = matrix(1,length(ly),3)
               X[,2] = t(lx[1,])-out1[j]
               X[,3] = t(lx[2,])-out2[i]
               if(!is.null(count)){
                   temp = temp*count[ind]
                   W1 = diag(temp, length(temp), length(temp))
               }else{
                   W1 = W
               }
              
               beta = ginv(t(X)%*%W1%*%X)%*%t(X)%*%W%*%ly
               rm(X,W,W1)
               mu[i,j] = beta[1]
           }else{
               invalid = 1
               cat("Not enough points in local window, please increase bandwidth!\n")
               return(list(invalid = invalid, mu = NULL))
           }
        }
    }
    
    if(!is.null(mu)){
       #mu[lower.tri(mu)] = mu[upper.tri(mu)] #assign lower triangular part of the mu matrix
                                             #to be the same as the upper triangular part
       a = matrix(0, dim(mu)[1], dim(mu)[2]);
       a[upper.tri(a)] = mu[upper.tri(mu)]
       mu = diag(mu)*diag(1, dim(mu)[1],dim(mu)[2])+a+t(a)
    }
    return(list(invalid = invalid, mu = mu))

}
