getVal = function(X, varname){

    id = match(varname, names(X));
    if(is.na(id)){
       Xname = as.character(substitute(X))
       cat(paste('Error: "', varname, '" is an invalid variable name! Use function names(', Xname, ') to get a valid name!\n', sep = ""));
       val = NULL
    }else{
       val = X[[id]]
    }
    val
}
