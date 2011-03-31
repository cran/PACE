interp11 = function(xx,y,newx = xx, method = "natural", nder = 0){
#xx is a vector, where all rows of y are evaluated at
#y is a matrix
   apply(y,2,function(x) interp1(xx,x,newx, method = method, nder = nder))

}