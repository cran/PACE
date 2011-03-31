mu_true = function(tt,p){
   
   if(length(p) == 1){
      lb = 0;
      ub = p;
   }else{
      lb = p[1];
      ub = p[2];
   }

   tt[!(tt >= lb & tt <= ub)] = 0
   mu = tt+sin(tt)
   mu

}