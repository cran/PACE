FPCAeval = function(yy, subjectID = NULL, newt){

   n = dim(getVal(yy,"y"))[1]
   if(is.null(subjectID)){
      subjectID = 1:n
   }
   if(any(subjectID > n | subjectID < 1)){
      stop("subjectID must be a vector of integer(s) that is between 1 and ", n)
   }
   
   if(is.matrix(newt)){
      sameTIME = 0;
      if(length(subjectID) != dim(newt)[1]){
         cat("Error: length of the subjectID and the row dimension of newt must be the same!\n")
         return(NULL)
      }
      
   }else if(is.vector(newt)){
      sameTIME = 1;
      newt = newt[!is.na(newt)]
   }else{
      stop("newt must be either a matrix or a vector!")
   }

   out1 = getVal(yy,"out1copy")
   no_opt = getVal(yy, "no_opt")
   mu = getVal(yy, "mucopy")
   phi = getVal(yy,"phicopy")
   if(is.matrix(phi)){
      phi = phi[,1:no_opt]
   }
   xi_est = getVal(yy,"xi_est")


   if(sameTIME){
     newmu = interp1(out1, mu, newt)
     newphi = interp11(out1, phi, newt);
     if(no_opt == 1){
        if(!is.matrix(newphi))
           phi = matrix(newphi, length(newt), no_opt);
     }
     MU = matrix(rep(newmu, length(subjectID)), length(subjectID),byrow = TRUE)
     y_pred = MU+xi_est[subjectID,]%*%t(newphi)

   }else{

        ni = sapply(1:length(subjectID), function(x) sum(!is.na(newt[x,])))
        y_pred = matrix(NA, length(subjectID), max(ni))

        for(i in 1:length(subjectID)){
           ti = newt[i,]
           idx = which(!is.na(ti))
           ti = ti[idx]
           mu_i = interp1(out1, mu, ti)
           phii = interp11(out1, phi, ti);
           if(no_opt == 1){
              if(!is.matrix(phii))
                 phii = matrix(phii, length(ti), no_opt);
           }
           y_pred[i,idx] = mu_i+xi_est[subjectID[i],]%*%t(phii)
        }
   }

   y_pred
 
}