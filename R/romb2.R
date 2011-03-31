romb2 = function(x,y,decdigs = 10){
#perform romberg integration for each column of y
  if(is.matrix(y)){
     m = dim(y)[2]
  }else{
     res = romb(x,y,decdigs)
     return(res)
  }
  res = numeric(m)
  for(i in 1:m){
    res[i] = romb(x,y[,i],decdigs)
  }
  res
}