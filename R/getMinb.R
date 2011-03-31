getMinb = function(tt,out1,regular = 0, npoly = 1){

   if(regular == 0){

      res = designPlotCount(tt,out1,1,1);   #obtain the count matrix based on the design plot 
      dstar = minb(out1,2+npoly);           #rough initial estimate of dstar based on 1-D sense
                                            #for at least 3 points
      res[,c(1,ncol(res))] = 1;
      res = t(res);
      ids = (res > 0);
      rm(res, tt)
      b = matrix(rep(out1,length(out1)),,length(out1));
      dstar = max(dstar,max(diff(b[ids]))/2);
      rm(b)

   }else if(regular == 1){

      dstar = minb(out1,1+npoly)*2;

   }else{

      dstar = minb(out1,2+npoly)*1.5;

   }

   dstar
}