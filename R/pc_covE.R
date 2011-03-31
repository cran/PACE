pc_covE = function(tt,out1,bw_xcov,ngrid = 51,cut = 1,kernel = c("gauss","epan","rect", "quar"),rcov, npoly = 1){

    kernel = kernel[1]
    ttt = vecData(tt)
    a0 = min(out1)
    b0 = max(out1)
    lint = b0 - a0
    h = (b0-a0)/(ngrid-1)
    out21 = seq(a0,b0,len = ngrid)
    out22 = out21

    tpairn = rcov$tpairn
    #tneq = apply(tpairn,2,function(x) x[1] != x[2])
    tneq = tpairn[1,] != tpairn[2,];
    cyy = rcov$cyy

    #This is for the case when regular = 1, the raw covariance
    #matrix needs to be divided by the number of individual sums
    #for each element in the matrix. In the case of regular = 2,
    #the division is n for each of the element.
    if(!is.null(rcov$count))
       cyy = cyy/rcov$count

    cxx = cyy[tneq]
    rm(rcov)

    win1 = rep(1,len = length(cxx))
    #get smoothed variance function for y(t) using lwls
    teq = !tneq
    vyy = cyy[teq]
    win2 = rep(1,len = length(vyy))
    #yvar is the variance function
   
    r = lwls(bw_xcov[1], kernel, 1, npoly, 0, tpairn[1,teq], vyy, win2, out21, 0)
    yvar = r$mu
    invalid = r$invalid

    if(invalid == 0){
       #estimate variance of measurement error term
       #use quadratic form on diagonal to estimate Var(x(t))
       r = rotate_mlwls(bw_xcov, kernel, tpairn[,tneq], cxx, win1, rbind(out21,out22),npoly)
       invalid = r$invalid
       xvar = r$mu
       rm(r)

       if(invalid == 0){
           if(cut == 0){
              sigma = romb(out21,yvar-xvar)/lint
           }else if(cut == 1){
              a = a0+lint*0.25
              b = a0+lint*0.75
              ind1 = out21 > a & out21 < b
              yvar1 = yvar[ind1]
              xvar1 = xvar[ind1]
              sigma = romb(out21[ind1],yvar1-xvar1)*2/lint      
           }
       }

    }
    if(sigma < 0){
       cat("Estimated sigma is negative, reset to zero now!\n")
       sigma = 0
    }

    list(invalid = invalid, sigma = sigma, xvar = xvar, yvar = yvar)

}
