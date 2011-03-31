getBins = function(x,y,xx,N){

    count = numeric(N-1)
    newy = count
   
    for(i in 2:(N-1)){
       ids = (x >=xx[i-1] & x<xx[i]);
       if(all(ids == 0)){
          count[i-1] = 0;
          newy[i-1] = NA;
       }else{
          count[i-1] = sum(ids);
          newy[i-1] = mean(y[ids])
       }
    }
    
    #for the last bin, it includes left and right end point
    ids = (x >=xx[i] & x <= xx[i+1])
    if(all(ids == 0)){
       count[i] = 0;
       newy[i] = NA;
    }else{
       count[i] = sum(ids);
       newy[i] = mean(y[ids]);
    }
    return(data.frame(newy = newy, count = count))
}
