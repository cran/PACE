vecData = function(x){
#row-wise vectorization, meaning vectorize data from each subject
#remove any NA from each subject
   x = as.matrix(x)
   x = t(x)
   x[!is.na(x)]

}