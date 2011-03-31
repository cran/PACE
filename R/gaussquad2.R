gaussquad2 = function(x,y,n = 30){

  if(is.matrix(y)){
     m = dim(y)[2]
  }else{
     res = gaussquad(x,y,n)
     return(res)
  }
  res = numeric(m)
  for(i in 1:m){
    res[i] = gaussquad(x,y[,i],n)
  }
  res


}