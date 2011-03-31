getTimeID = function(tt, n, isRandom){

   ni = numeric(n)
   tjID = ni

   if(isRandom == 1){
       for(i in 1:n){
          ni[i] = sum(!is.na(tt[i,]))
          if(ni[i] > 1)
             tjID[i] = sample(1:ni[i],1)     #find random index of random time point to be used in the CV prediction par
                                             #for ni >= 2
       }
   }else{
       for(i in 1:n){
           ni[i] = sum(!is.na(tt[i,]))
           tjID[i] = ((1000+i) %% ni[i])+1  
       }
   }
   
   list(tjID = tjID, ni = ni)

}