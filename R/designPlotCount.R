designPlotCount = function(tt,out1,noDiagonal = 1, isIndicator = 1){

   N = length(out1) #number of unique time points
   res = matrix(0, N,N)

   for(i in 1:dim(tt)[1]){
      cur = tt[i,]   #time points from subject i
      cur = cur[!is.na(cur)]
      curidx = which(out1 %in% cur);
      if(isIndicator){
          res[curidx,curidx] = 1
      }else{

          res[curidx,curidx] = res[curidx,curidx]+1
      }
   }   
   if(noDiagonal)
      diag(res) =0
   res

}