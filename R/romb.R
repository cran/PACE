romb = function(x,y,decdigs = 10){

  rom = matrix(0,2,decdigs)
  romall = numeric(2^(decdigs-1)+1)
  a = min(x)
  b = max(x)
  romall = interp1(x,y,seq(a,b,by = (b-a)/2^(decdigs-1)))
  h = b-a
  rom[1,1] = h*(romall[1]+romall[length(romall)])/2
  for(i in 2:decdigs){
      st=2^(decdigs-i+1);
      rom[2,1]=(rom[1,1]+h*sum(romall[st/2+seq(1,2^(decdigs-1), by = st)]))/2
    
      for(k in 1:(i-1)){
         rom[2,k+1]=((4^k)*rom[2,k]-rom[1,k])/((4^k)-1)
      }
      rom[1,1:i]=rom[2,1:i]
      h=h/2
  }
  res = rom[1,decdigs]
  return(res)
}
