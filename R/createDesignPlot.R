createDesignPlot = function(datafile, isColorPlot = 0, noDiagonal= 1, isIndicator = 1, yname = "y"){

    if(is.character(datafile)){
      timeMatrix = read.table(datafile, fill = TRUE) #time matrix, where missing is denoted by NA
      #tt = convertData(timeMatrix)      #convert it to a list
    }else{ 
      tt = datafile                      #no conversion is needed, if the input is already a matrix or data.frame
    }
    
    out1 = unique(sort(tt[!is.na(tt)]))
    res = designPlotCount(tt,out1, noDiagonal, isIndicator)
    
    createBlackPlot = function(res, out1, yname = "y"){
      x11();
      plot(0,xlab = expression(T[ij]), ylab = expression(T[ik]),xlim = c(min(out1),max(out1)),
           ylim = c(min(out1),max(out1)),main = paste("Design Plot of", yname), font = 2, type = "n")  
      for(i in 1:length(out1)){
         idx = which(res[i,] > 0)
         points(rep(out1[i],len = length(idx)),out1[idx],type = "p",pch = 20, col = "black")
      }   
    
      invisible(NULL)

    }

    createColorPlot = function(res,out1,yname = "y"){

       x11()

       searchCol = function(val){

          if(val == 1){
              col = "red"
          }else if(val == 2){
              col = "magenta"
          }else if(al >=3 && val <= 5){
              col = "green"
          }else if(val > 6){
              col = "blue"
          }
          col
       }

      plot(0,xlab = expression(T[ij]), ylab = expression(T[ik]),xlim = c(min(out1),max(out1)),
           ylim = c(min(out1),max(out1)),main = paste("Design Plot of", yname), font = 2, type = "n")  
      for(i in 1:length(out1)){
         tmp = res[i,]
         idx = which(tmp > 0)
         for(j in 1:length(idx)){
            points(out1[i],out1[idx[j]],type = "p",pch = 20, col = searchCol(tmp[idx[j]]))
         }
      }   
      invisible(NULL)
    }

    if(isColorPlot){
       createColorPlot(res,out1,yname)
    }else{
       createBlackPlot(res,out1,yname)
    }


}
