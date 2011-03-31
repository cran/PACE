
FPCApred = function(yy, newy, newt, regular = NULL){

   p = getVal(yy,"ops");
   if(is.null(regular)){
      regular = getVal(p, "regular");
   }
   
   nsub = dim(newy)[1]

   if(is.vector(newt)){
      newt = newt[!is.na(newt)];
      newt = matrix(rep(newt, nsub), nsub,byrow = TRUE)
      regular = 2;
   }

   mu = getVal(yy,"mucopy");
   no_opt = getVal(yy,"no_opt");
   phi = getVal(yy,"phicopy");
   if(no_opt == 1){
      if(!is.matrix(phi))
           phi = matrix(phi, length(phi), no_opt);
   }else if(no_opt > 1){
      phi = phi[,1:no_opt];
   }
   lambda = getVal(yy,"lambda")
   lambda = lambda[1:no_opt];
   sigma = getVal(yy,"sigma");
   sigmanew = getVal(yy,"sigmanew")
   
   error = getVal(p, "error");
   method = getVal(p, "method");
   shrink = getVal(p, "shrink");
   out1 = getVal(yy, "out1copy");
   rho = getVal(yy, "rho_opt");     

   #rr$ypred, rr$xi_new and rr$xi_var 
   rr = getXI(newy, newt, mu, phi, lambda, sigma, sigmanew, no_opt, error, method, shrink, out1, regular, rho);
   rr

}
