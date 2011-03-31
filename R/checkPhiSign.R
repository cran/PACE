checkPhiSign = function(out1, phi, true_phi){

  a1 = romb2(out1, (phi-true_phi)^2)
  a2 = romb2(out1, (-phi-true_phi)^2)

  if(a1 <= a2){
    return(1)
  }else{
    return(-1)
  }
  
}
