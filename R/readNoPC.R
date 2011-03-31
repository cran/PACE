readNoPC = function(FVE, verbose = "on"){
       
       cat(paste("Enter the number of principal components you want to choose (K <= ", length(FVE),"):\nK=", sep = ""));
       no_opt = as.numeric(readLines(n = 1))

       if(no_opt > length(FVE)){
             cat(paste("K can't be larger than ", length(FVE), "! Reset to ", length(FVE), " now!\n", sep = ""))
             no_opt = length(FVE)
       }
       if(verbose == "on"){
          cat(paste("You just chose ", no_opt,  " principal component(s).\n", sep = ""));
          cat(paste("It accounts for ", round(FVE[no_opt], digits = 4)*100, "% of total variation.\n", sep = ""));
       }

       invisible(no_opt)
}