mapX1d = function(xx,yy,newx){

   id = sapply(newx, function(x) which(x == xx))
   if(is.vector(yy)){
        yy[id]
   }else if(is.matrix(yy)){
        yy[id,]
   }else{
        stop("y cannot be empty!");
   }
}