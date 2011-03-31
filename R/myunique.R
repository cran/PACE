myunique = function(xx, sorted = TRUE){
#get the unique value or unique sorted value of input vector "xx"
#any missing values from "xx" are ignored
#out1[id] is same as "xx" (without missing values)

   xx = xx[!is.na(xx)]
   if(sorted){
     out1 = unique(sort(xx))
   }else{
     out1 = unique(xx)
   }
   id = sapply(xx, function(x) which(out1 == x))
   list(out1 = out1, id = id)

}