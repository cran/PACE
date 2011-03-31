getEigens = function(xcov, out1, out21, noeig, disp = FALSE){
  r = range(out21)
  h = (r[2]-r[1])/(length(out21)-1)
  
  ngrid = dim(xcov)[1]
    
  r = eigen(xcov)
  eigens = r$vectors
  d = r$values
  rm(r)
  
  idx = which(Im(d) != 0)      #find indices for imaginary eigenvalues
  if(length(idx) >  0){
      stop(paste(length(idx), "eigenvalues are complex. The estimated auto-covariance surface is not symmetric!"));
  }
  idx = which(d <= 0)
  if(length(idx) > 0){
       if(disp)
          cat(paste(length(idx), "real eigenvalues are negative or zero and are removed!\n"))
       eigens = eigens[, d > 0]
       d = d[d > 0]
  }

 
  if(noeig > length(d)){
    noeig = length(d)
    cat(paste("At most", noeig, "number of PC can be selected!\n"))
  }

  eigens = eigens[,1:noeig]
  if(!is.matrix(eigens))
     eigens = matrix(eigens,length(out21),noeig)
  d = d[1:noeig]
  eigens = eigens/sqrt(h)
  lambda = h*d
  
  for(i in 1:noeig){
      eigens[,i] = eigens[,i]/sqrt(romb2(out21,eigens[,i]^2))
      if(eigens[2,i] < eigens[1,i])
        eigens[,i] = -eigens[,i]
  }
  
  #interpolate from the normalized the eigenfunctions
  phi = interp11(out21,eigens, out1)

  #normalize smoothed eigenfunctions
  for(i in 1:noeig){
      phi[,i] = phi[,i]/sqrt(romb2(out1,phi[,i]^2))
  }

  list(lambda = lambda, phi = phi, eigen = eigens, noeig = noeig)
  
}
