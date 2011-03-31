minb = function(x, numPoints = 2){

   n = length(x)
   x = sort(x)
   if(numPoints > 1){
     max(x[numPoints:n]-x[1:(n-numPoints+1)])  
   }else{
     max((x[2:n]-x[1:(n-1)])/2)
   }

}
