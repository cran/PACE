setVal = function(X, varname, val){

   id = match(varname, names(X));
   if(is.na(id)){
       Xname = as.character(substitute(X))
       cat(paste('Error: "', varname, '" is an invalid variable name! Use function names(', Xname, ') to get a valid name!\n', sep = ""));
   }else{
       X[[id]] = val;
   }
   return(X)
}