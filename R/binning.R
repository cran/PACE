binning = function(x,y,M = 10, isMnumBin = 1, nonEmptyOnly = 0, a0=min(x, na.rm = TRUE), b0 = max(x,na.rm = TRUE)){


    if(M <= 0){
      stop("M must be positive!")
      invisible(NULL)
    }
    
    x = x[!is.na(x)]
    y = y[!is.na(y)]

    if(isMnumBin == 0){
       h = M;
       rx = range(x, na.rm = TRUE)
       rx = rx[2]-rx[1]
       if(h>=rx){
           res = getResMisOne(x,y,h);
           return(res)
       }
       xx = seq(a0,b0,by = h)
       M = length(xx)-1;
 
    }else if(isMnumBin){       
       if(M == 1){
          res = getResMisOne(x,y);
          return(res)
       }
       xx = seq(a0,b0, len = M+1);
       h = xx[2]-xx[1]
    }
    
    N = length(xx);
    midpoint = xx[1:(N-1)]+h/2;
 
    #xx[i-1] is like the lower limit
    #xx[i] is like the upper limit
    #in each bin, it includes the left end point
    
    r = getBins(x,y,xx,N);
    newy = r$newy
    count = r$count

    if(nonEmptyOnly){

        midpoint =midpoint[count > 0]
        newy = newy[count>0]
        count = count[count>0]
        M  = length(midpoint)

    } 
    
    res = list(midpoint = midpoint, newy = newy, count = count, 
               numBin = M, binWidth = h)
    return(res)

}
