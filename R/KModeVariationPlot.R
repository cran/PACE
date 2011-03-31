KModeVariationPlot = function(yy, k = 1:getVal(yy,"no_opt")){

  no_opt = getVal(yy,"no_opt");
  if(any(k <= 0)){
     cat("Error: k must be positive integer(s)!\n");
     return(NULL);
  }else if(any(k > no_opt)){
     cat(paste("Error: any k cannot be larger than no_opt = ", no_opt, "\n", sep = ""));
     return(NULL);
  }  

  lambda = getVal(yy,"lambda");
  mu = getVal(yy,"mu");

  if(is.list(mu)){
    mu = mu[[1]];
  }  
  
  out1 = getVal(yy,"out1");
  phi = getVal(yy,"phi");
  
  if(is.list(phi)){
     phi = phi[[1]];
  }

  for(i in 1:length(k)){
     x11();
     alpha = c(-2,-1,0,1,2)*sqrt(lambda[k[i]]);
     cmat =  rainbow(length(alpha),start = 0.5, end = .9, alpha = .5);
     tmp = matrix(NA, 5, length(mu));
     for(j in 1:length(alpha)){
        tmp[j,] = mu+alpha[j]*phi[,k[i]];       
     }
     minY = min(tmp);
     maxY = max(tmp);
     plot(out1, tmp[1,], xlab = "t", ylab = parse(text = paste("mu(t)+alpha*phi[",k[i], "](t)", sep = "")), 
           main = paste("K-th mode of variation plot (k = ", k[i], ")",  sep = ""), type = "l", col = cmat[1], xlim = c(min(out1),max(out1)),
           ylim = c(minY, maxY));
     for(j in 2:length(alpha)){
        lines(out1,tmp[j,], col = cmat[j]);
     }
     legend(x = "topleft", legend = c(expression(alpha == -2),expression(alpha == -1),
            expression(alpha== 0), expression(alpha == 1), expression(alpha == 2)), col = cmat, lty = 1, inset = .01, box.col = "white", cex = 0.8)
  } 

}
