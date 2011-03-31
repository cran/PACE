interp1 = function(x,y,newx = x, method = "natural", nder = 0){
#use splinefun() in the stats package of R to perform 1-D interpolation
#default is used natual cubic spline, other possible values, see
#the description for splinefun().    
    f = splinefun(x,y,method = method)
    if(is.matrix(newx)){
      res = f(newx[,1], deriv = nder) 
      cols = ncol(newx)
      for(i in 2:cols){
        res = cbind(res, f(newx[,i],deriv = nder))
      }
   }else{
     res = f(newx,deriv = nder)
   }
   res

}