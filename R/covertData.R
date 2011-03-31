#convert the input data, possibly had missing values
#in the form of data.frame or matrix to a list in R
convertData = function(x){
   xx = apply(x, 1,function(x) list(x))
   x = lapply(xx,function(x) {
                   tmp = x[[1]]; 
                   tmp[!is.na(tmp)]}
             )      
   x 
}