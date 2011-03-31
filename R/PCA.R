library(akima) #for interpolations

ask = function(){
  cat("Press any key to continue ...")
  invisible(readLines(n = 1))
}

PCA = function(y,tt, bwmu = 0, bwmu_gcv = 1, bwxcov = c(0,0) ,bwxcov_gcv = 1, 
               ntest1 = 30, ngrid1 = 30, selection_k = "BIC1", FVE_threshold = 0.85, 
               maxk = 20, control = "auto", regular = 0, error = 1, ngrid = 51, 
               method = "CE", shrink = 0, newdata = NULL, kernel = "gauss", 
               numBins = NULL, yname = NULL, screePlot = 0, designPlot = 0, rho = "cv", 
               verbose = "on")
{

   if(is.data.frame(y)){
      y = as.matrix(y)
   }
   if(is.data.frame(tt)){
      tt = as.matrix(tt)
   }
   cut = NULL;  #use default value of cut
   if(is.null(yname)){
     yname = as.character(substitute(y))
   }

   if(is.null(bwmu)){
      bwmu = 0;   #bandwidth choice for mean function is using CV or GCV
   }

   if(is.null(bwmu_gcv)){
      bwmu_gcv = 1;  #bandwidth choice for mean function is GCV if bwmu = 0
   }

   if(is.null(bwxcov)){
      bwxcov = c(0,0); #bandwidth choice for covariance function is CV or GCV
   }

   if(is.null(bwxcov_gcv)){
      bwxcov_gcv = 1;  #bandwidth choice for covariance function is GCV if bwxcov = c(0,0)
   }

   if(is.null(ntest1)){
      ntest1 = 30;
   }

   if(is.null(ngrid1)){
      ngrid1 = 30;
   }

   if(is.null(selection_k)){
      selection_k = "BIC1";
   }

   if(selection_k == "FVE" && is.null(FVE_threshold)){
      FVE_threshold = 0.85;
   }

   if(is.null(error)){
      error = 1;    #error assumption with measurement error
   }

  if(is.null(maxk)){
     maxk = 20;    #20 PC candidates
  }

  if(is.null(shrink)){
    shrink = 0;
  }

  if(is.null(method)){
     shrink = 0;
     method = "CE";   #method to estimate the PC score is through conditional expectation
  }

  if(is.null(regular)){
     regular = 0;     #sparse functional data
  }

  if(is.null(kernel)){
     if(regular == 2){
        kernel = "epan";   #kernel: Epanechnikov
     }else{
        kernel = "gauss";  #kernel: Gaussian
     }
  }
  
  if(shrink == 1 && (error != 1 || method != "IN")){
     cat('shrinkage method only had effects when method = "IN" and error = 1! Reset to shrink = 0 now!\n');
     shrink = 0      
  }

  if(is.null(control)){
     control = "auto"; #automatically choose the number of PC not through interactive input
                       #by the user after visualization
  }
  
  if(is.null(cut) && error == 1){
      cut = 1;         #cut off the boundary
  }else if(is.null(cut)){
      cut = 0;
  }

  if(error == 0 && cut == 1){
     cat("cut is only effective if error = 1! Reset to cut = 0 now!\n");
     cut = 0;
  }

  if(is.null(ngrid)){
    ngrid = 51;    #number of output time grids is 51
  }

  if(maxk > (ngrid-2)){
    cat(paste("maxk can only be less than or equal to", ngrid-2,"! Reset to be", ngrid-2, "now!\n"));
    maxk = ngrid -2;
  }

  if(is.numeric(selection_k)){

    if(selection_k > (ngrid-2)){
        cat(paste("maxk can only be less than or equal to", ngrid-2,"! Reset to be", ngrid-2, "now!\n"));
        maxk = ngrid -2;
    }else if(selection_k <= 0){
       cat("selection_k must be a positive integer! Reset to BIC1 now!\n");
       selection_k = "BIC1"
       FVE_threshold = 0.85;
    }
  }

  if(is.null(screePlot)){
     screePlot = 0;
  }

  if(is.null(designPlot)){
     designPlot = 0;
  }

  if(is.null(rho)){
     rho = "cv";
  }
  
  if(is.null(verbose)){
     verbose = "on";
  }

  if(error == 0 && (selection_k == "AIC2" || selection_k == "BIC2")){
     cat('When assume no measurement error, cannot use "AIC2" or "BIC2". Reset to "BIC1" now!\n')
     selection_k = "BIC1"
  }

  kernNames = c("rect", "gauss", "epan", "gausvar", "quar");
  if(!(kernel %in% kernNames)){
     cat(paste('kernel', kernel, 'is unrecognizable! Reset to "gauss" now!\n'));
     if(regular == 2){
         kernel = "epan";
     }else{
         kernel = "gauss";
     }
  }  

  ncohort = nrow(tt);  #obtain the number of curves or subjects
  ni = apply(tt,1,function(x) sum(!is.na(x)))

  if(all(ni == 1)){
     cat("Error:FPCA is aborted because the data do not contain repeated measurements!\n"); 
     return(NULL);    
  }   

  #Prebin the data if conditions are satisfied as specified in the help for numBins
   
  if(is.null(numBins)){
     res = binData(y, tt, regular, verbose);
     newy = res$newy;
     newt = res$newt;
     rm(res);
  }else if(numBins >= 10){
     res = binData(y,tt, regular, verbose, numBins);
     newy = res$newy;
     newt = res$newt;
     rm(res);
  }else if(numBins == 0){
     newy = NULL;     #no binning set by user
  }else if(numBins < 10 || numBins < 0){
     newy = NULL      #no binning due to number of bins is too small
     cat("number of bins must be at least 10! No binning will be performed!\n");
  }

  if(!is.null(newy)){
     y = newy;
     tt = newt;
     if(regular == 0){
        regular = 1;
     }
     rm(newy, newt)
  }

  ops = list(bwmu = bwmu, bwmu_gcv = bwmu_gcv, bwxcov = bwxcov, bwxcov_gcv =  bwxcov_gcv, ntest1 = ntest1,
             ngrid1 = ngrid1, selection_k = selection_k,  FVE_threshold = FVE_threshold, maxk =  maxk,
             control = control, regular = regular, error = error, ngrid = ngrid, method = method, shrink = shrink,
             newdata = newdata, kernel = kernel, numBins = numBins, yname = yname, screePlot = screePlot,
             designPlot = designPlot, rho = rho, verbose = verbose);

  #pool all the subjects and their corresponding time points into 1 x N vectors
  #ttt: vector to hold time points
  #yy: vector to hold the observed measurements
  ttt = vecData(tt)    # 1 x N vector to hold the observed time points from all subjects
   yy = vecData(y);    # 1 x N vector to hold the observed measurements from all subjects

  #Initial out1 is based on the unique time points of the pooled data + the unique
  #time points of "newdata", the output time grid. When newdata = [], output
  #"out1" is equivalent to be the unique sorted pooled time points; otherwise, it
  #corresponds to the unique "newdata".
  out1 = myunique(c(ttt, newdata))$out1;

  if(designPlot == 1){
      createDesignPlot(tt, 0, 1, 1, yname)    
  }

  if(verbose == "on"){
     cat("Part I: Obtain smoothed mean curve\n");
  }

  #when bwmu = 0 and bwmu_gcv = 0, use leave-one-curve-out CV method for bw choice
  #when bwmu = 0 and bwmu_gcv = 1, use leave-one-out GCV method for bw choice
  #when bwmu > 0, user-defined bw for mean function
   
  if(bwmu == 0){
      if(bwmu_gcv == 1){

         bw_mu = gcv_lwls(yy,ttt, kernel, 1,1,0,regular,verbose)$bopt;  #use GCV method to choose bw for mean function
         if(is.na(bw_mu)){
            return(NULL);
         }

         bw_mu = adjustBW1(kernel, bw_mu, 1,0, regular, verbose);

      }else{                                                            #use CV method to choose bw for mean function
         bw_mu = cvfda_lwls(y,tt,kernel,1,1,0,regular,verbose);
      }
  }else if(bwmu > 0){
      bw_mu = bwmu;
  }else{
      cat("Error: Bandwidth choice for the mean function must be positive!\n");
      return(NULL);
  }

  #define the vector of case weight in the local weighted least square
  #here, it is set to be one for all subjects
  win1 = rep(1, len = length(ttt))

  rr = lwls(bw_mu, kernel, 1, 1, 0, ttt, yy, win1, out1)
  invalid = rr$invalid;
  mu = rr$mu;
  rm(rr)

  if(verbose == "on"){
     cat("Part II: Choose bandwidth of smoothing covariance surface\n");
  }

  rcov = getRawCov(y,tt,out1, mu, regular, 0);   #obtain raw covariance
  
  if(bwxcov[1] == 0 || bwxcov[2] == 0){
      if(bwxcov_gcv == 1){
          bw_xcov = gcv_mullwlsn(tt, ngrid1, regular, error, kernel, rcov, verbose)$bw_xcov;
          
          if(any(is.na(bw_xcov))){
             cat("Error: FPCA is aborted because the observed data is too sparse to estimate the covariance function!\n");
             return(NULL);
          }
          bw_xcov = adjustBW2(kernel, bw_xcov, 1, 0, regular, verbose);

      }else{
          bw_xcov = cv_mullwlsn(y,tt,mu,ntest1,ngrid1,regular, error, kernel,rcov, verbose);
      }
  }else if(all(bwxcov > 0)){
      bw_xcov = bwxcov; 
  }else if(bwxcov[1] < 0 || bwxcov[2] < 0){
      cat("Error: Bandwidth choice for the covariance function must be positive!\n");
      return(NULL);
  }

  if(verbose == "on"){
     cat("Part III: Choose number of principal components functions\n");
  }

  AB_method = c("full", "rand");
  out21 = seq(min(out1), max(out1), len = ngrid);
   
  rcov1 = rcov;

  if(error == 1){
     tpairn = rcov1$tpairn;
     tneq = tpairn[1,] != tpairn[2,];
     cyy = rcov1$cyy;
     rcov1$tpairn = tpairn[,tneq]
     rcov1$cxxn = cyy[tneq];
     rcov1$win = rep(1,len = length(rcov1$cxxn));
     if(regular == 1){
        rcov1$count = rcov1$count[tneq];
     }
  }

  if(regular == 1){  #smooth raw covariance
     rr = mullwlsk(bw_xcov,kernel,rcov1$tpairn,rcov1$cxxn,rcov1$win,out21,out21,rcov1$count) 
  }else{             #smooth raw covariance
     rr = mullwlsk(bw_xcov,kernel,rcov1$tpairn,rcov1$cxxn,rcov1$win,out21,out21)
  }
  invalid = rr$invalid
  xcov = rr$mu
  rm(rr)

  xcov = (xcov+t(xcov))/2; #transform the smoothed covariance matrix to guarantee it is a symmetric matrix.

  if(invalid == 0){
     rr = no_FVE(xcov, FVE_threshold);     
     no_opt = rr$no_opt;
     FVE = rr$FVE;
     rm(rr);
  }else{
     cat("FPCA is aborted because enough points to estimate the smooth covariance function!\n")
     return(NULL)
  }
  
  no_optCopy = no_opt;

  pc_options = c("AIC1", "AIC2", "BIC1", "BIC2", "FVE", "user", "AIC_R");
  AIC = NULL;
  BIC = NULL;

  if(is.character(selection_k)){
     k_id = match(selection_k,pc_options);
     if(is.na(k_id)){
        cat(paste('Invalid method name for selection_k! Reset to "FVE" method with threshold =', FVE_threshold,'!\n'));
        k_id = 5;
     }
     if(k_id == 1 || k_id == 2){
        rr = no_AIC(y,tt,mu, bw_xcov, ngrid, regular, maxk, AB_method[k_id],
                 method, shrink, out1, out21, kernel, error, cut, rcov, xcov)
        no_opt = rr$no_opt;
        AIC = rr$aic;
        rm(rr);        
     }else if(k_id == 3 || k_id == 4){
        rr = no_BIC(y,tt, mu, bw_xcov, ngrid, regular, maxk, AB_method[k_id-2],
                  method, shrink, out1, out21, kernel, error, cut, rcov, xcov);
        no_opt = rr$no_opt;
        BIC = rr$bic;
        rm(rr)
     }else if(k_id == 7){
        no_opt = ngrid;
     }
     if(is.null(no_opt)){
        k_id = 5;
        no_opt = no_optCopy;
     }     

  }else if(is.numeric(selection_k) && selection_k > 0){
     no_opt = selection_k;
     k_id = 6;
  }else{
     cat(paste('"selection_k" must be a positive integer! Reset to "FVE" method with threshold =', FVE_threshold, '!\n'));
     k_id = 5; 
  }

  if(k_id != 7){
     if(verbose == "on"){
        cat(paste("Best number of principal components selected by ", pc_options[k_id]," : ", no_opt, ".\n", sep = ""));
     }
     if(k_id != 5){
        if(verbose == "on"){
           cat(paste("It accounts for ", round(FVE[no_opt], digits = 4)*100, "% of total variation.\n", sep = ""));
        }
     }else{
        if(verbose == "on"){
           cat(paste("It accounts for ", round(FVE[no_opt], digits = 4)*100, "% of total variation (threshold = ", FVE_threshold, ").\n", sep = ""));
        }
     }

     if(no_opt == maxk && k_id < 5){
        cat(paste(pc_options[k_id], " cannot find the best No. of PC for maxk = ", maxk, ". Increase maxk to get better results.\n", sep = ""))
     }

     if(verbose == "on"){
        cat(paste("FVE calculated from ", ngrid, " possible eigenvalues: \n", sep = ""))
        print(FVE)
     }

     if(control == "look"){
        no_opt = no_SP(FVE, ngrid, no_optCopy, FVE_threshold, yname, verbose);
     }

     #Now output scree plot based on final no_opt
     if(screePlot == 1){
         createSP(FVE, no_opt, yname)         
     }

  }
  
  if(verbose == "on"){
     cat("Part IV: Perform principal components analysis\n");
  }

  if(error == 1){
     rr =pc_covE(tt,out1,bw_xcov,ngrid,cut,kernel,rcov);
     invalid = rr$invalid;
     sigma = rr$sigma;
     rm(rr);

     if(invalid == 0){
         rr = pc_est(y, tt, mu, xcov, sigma, no_opt, error, method, shrink, out1, out21,regular,rho, verbose);
         xi_est = rr$xi_est;
         xi_var = rr$xi_var;
         lambda = rr$lambda;
         phi = rr$phi;
         eigens = rr$eigen;
         no_opt = rr$noeig;
         xcovfit = rr$xcovfit;
         y_predOrig = rr$y_predOrig;
         rho_opt = rr$rho_opt;
         sigmanew = rr$sig1;
         rm(rr);
     }else{
         xi_est = NULL; xi_var = NULL; lambda = NULL; phi = NULL;
     }
  }else if(error == 0){
     sigma = NULL;
     rr = pc_est(y, tt, mu, xcov, sigma, no_opt, error, method, shrink, out1, out21,regular,rho, verbose);
     xi_est = rr$xi_est;
     xi_var = rr$xi_var;
     lambda = rr$lambda;
     phi = rr$phi;
     eigens = rr$eigen;
     no_opt = rr$noeig;
     xcovfit = rr$xcovfit;
     y_predOrig = rr$y_predOrig;
     rho_opt = rr$rho_opt;
     sigmanew = rr$sig1;
     rm(rr);
  }

  if(verbose == "on"){
     cat("Part V: compute the smoothed individual trajectories\n");
  }

  #Save these copies for internal use in the FPCreg(), where the time points are guarantteed
  #to be evaluated at the distinct time points of pooled t's.
  mucopy = mu;
  phicopy = phi;
  eigenscopy = eigens;
  out21copy = out21;
  out1copy = out1;
  xcovcopy = xcov;
  xcovfitcopy = xcovfit;

  #Map the results to the user-defined output time grid "newdata"
  if(!is.null(newdata)){
     mu = interp1(out1, mu, newdata)
     out21 = seq(min(newdata), max(newdata), len = ngrid);
     phi = interp11(out1copy, phicopy, newdata)       #eigenfunction based on newdata
     eigens = interp11(out21copy, eigenscopy, out21)  #eigenfunction based on ngrid's of newdata

     if(no_opt == 1){
          if(!is.matrix(phi))
              phi = matrix(phi, length(out1), no_opt);
          if(!is.matrix(eigens))
              eigens = matrix(eigens,length(out21),no_opt)
     }
     o21 = expand.grid(out21copy, out21copy)
     xcovvec = as.vector(xcovcopy)
     o21new = expand.grid(out21,out21)
     xcov = interpp(o21[,1],o21[,2],xcovvec,o21new[,1], o21new[,2])$z
     xcov = matrix(xcov, length(out21), length(out21))
     xcovvec = as.vector(xcovfitcopy)
     xcovfit = interpp(o21[,1],o21[,2],xcovvec,o21new[,1],o21new[,2])$z
     xcovfit = matrix(xcovfit, length(out21), length(out21))
     out1 = newdata
  }

  if(invalid == 0){
     MU = matrix(rep(mu, dim(y)[1]), dim(y)[1],byrow = TRUE)
     y_pred = MU+xi_est%*%t(phi)
   }

  if(error == 0){
     sigma = NULL
  }

  list(no_opt = no_opt,sigma = sigma,lambda = lambda, phi = phi, eigens = eigens, xi_est = xi_est,
       xi_var = xi_var, mu = mu, bw_mu = bw_mu, xcov = xcov, bw_xcov = bw_xcov, xcovfit = xcovfit, 
       AIC = AIC, BIC = BIC,FVE = FVE, y_pred = y_pred, y_predOrig = y_predOrig, out1 = out1, out21 = out21,
       y = y, tt = tt, regular = regular, rho_opt = rho_opt, sigmanew = sigmanew, mucopy = mucopy, 
       phicopy = phicopy, eigenscopy = eigenscopy, out1copy = out1copy, out21copy = out21copy, 
       xcovcopy = xcovcopy, xcovfitcopy = xcovfitcopy, ops = ops)
  
}