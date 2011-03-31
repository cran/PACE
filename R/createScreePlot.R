createScreePlot = function(yy){

   ops = getVal(yy,"ops");
   yname = getVal(ops, "yname");
   createSP(getVal(yy,"FVE"), getVal(yy, "no_opt"), yname);
 
}