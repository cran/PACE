gaussquad = function(x,y,n = 30){

   a = min(x)
   b = max(x)
   
   res = lgwt(n,a,b)
   y_node = interp1(x,y,res$x)
   res = y_node%*%res$w
  
   if(dim(res)[1] == 1 && dim(res)[2] == 1)
     res = as.numeric(res)
   res

}
