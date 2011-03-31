getResMisOne = function(x,y,h = NULL){

    r = range(x, na.rm = TRUE)
    r = r[2]-r[1]
    if(is.null(h))
       h = r;
    M = 1
    midpoint = r/2
    count = length(x)
    newy = mean(y)
    res = list(midpoint = midpoint, newy = newy, count = count,
    numBin= M, binWidth = h)
    return(res)

}