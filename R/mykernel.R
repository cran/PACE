mykernel = function(x,kernel = "gauss"){
   if(kernel == "quar"){
      (15/16)*(abs(x) <= 1)*(1-x^2)^2
   }else if(kernel == "epan"){
      0.75*(abs(x)<=1)*(1-x^2)
   }else if(kernel == "rect"){
      0.5*(abs(x)<=1)
   }else if(kernel == "gausvar"){
      (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.25-0.25*x^2)
   }else if(kernel == "gausvar1"){
      (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.5-0.5*x^2)
   }else if(kernel == "gausvar2"){
      k1 = (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.25-0.25*x^2);
      k2 = (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.5-0.5*x^2);
      k1*k2
   }else{
      (1/sqrt(2*pi))*exp(-0.5*x^2);
   }

}
