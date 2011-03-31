xeig3 = function(tt, T){

  phi = matrix(0, 3, length(tt))
  if (all(tt >= 0 & tt <= T)){
     phi[1,] = -sqrt(2/T)*cos(2*pi*tt/T);
     phi[2,] = sqrt(2/T)*sin(2*pi*tt/T);
     phi[3,] = -sqrt(2/T)*cos(4*pi*tt/T);
  }else{
     cat("Error: t should be smaller than or equal to T!\n")
     return(NULL)
  }
  phi

}

