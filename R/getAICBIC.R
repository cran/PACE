getAICBIC = function(yy){

  p = getVal(yy,"ops");
  error = getVal(p, "error");
  maxk = getVal(p, "maxk");
  regular = getVal(p, "regular");
  method = getVal(p, "method");
  shrink = getVal(p, "shrink")

  y = getVal(yy, "y");
  tt = getVal(yy, "tt");

  aic1 =  rep(Inf, maxk)
  aic2 = aic1;
  bic1 = aic1;
  bic2 = aic1;

  N = length(vecData(tt))

  if(error == 1){
    sigma = getVal(yy,"sigma")
  }else{
    sigma = 0
  }
  
  out1 = getVal(yy, "out1copy");
  out21 = getVal(yy, "out21copy");
  xcov = getVal(yy, "xcovcopy");
  
  r = getEigens(xcov, out1, out21, maxk)
  lambda = r$lambda
  phi = r$phi
  rm(r)

  mu = getVal(yy, "mucopy");

  r1 = convertMuPhi(tt,out1, mu, phi, regular)
  muSub = r1$muSub
  phiSub = r1$phiSub
  rm(r1)

  for(k in 1:maxk){
  
       sig = lambda[1:k]
       
       r = getLogLik1(y,sig, muSub, phiSub, sigma, regular)
       invalid = r$invalid
       logLik = r$logLik
       rm(r)
       if(invalid){
             cat("The covariance matrix of the estimated function is nearly singular!\n")
       }else{
          aic1[k] = logLik+2*k;
          bic1[k] = logLik+log(N)*k;
       }

       logLik = getLogLik2(y,tt,sigma,k,method,shrink,regular,muSub, phiSub, sig)
       aic2[k] = logLik+2*k;
       bic2[k] = logLik + log(N)*k

  }
  
  list(aic1 = aic1, aic2 = aic2, bic1 = bic1, bic2 = bic2);

}
