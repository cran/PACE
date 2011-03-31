myrep = function(tt, each = TRUE){
#get a vector of repeated values
#ex: tt = matrix(c(1 ,2 ,3 ,NA,NA, seq(.1,.5,len = 5)),2,5, byrow = TRUE)
# myrep(tt): 1,1,1,2,2,2,3,3,3,.1,.1,.1,.1,.1,...,.5,.5,.5,.5,.5
# myrep(tt,FALSE): 1,2,3,1,2,3,1,2,3,.1,.2,...,.5,.1,...,.5,...,.1,...,.5
# where tt is a matrix or data.frame, possibly contains NA values
     if(each){
          xx1 = apply(tt,1, function(x) { 
                      tmp = x[!is.na(x)]
                      rep(tmp, each = length(tmp))
                      })
     }else{
          xx1 = apply(tt,1, function(x) { 
                      tmp = x[!is.na(x)]
                      rep(tmp, length(tmp))
                      })
     }
     as.numeric(unlist(xx1))

}
