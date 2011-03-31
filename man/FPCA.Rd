\name{FPCA}
\alias{FPCA}
\title{Principal Analysis by Conditional Expectation}
\description{A function to obtain functional principal components for sparsely or densely sampled 
             trajectory data without pre-smoothing of trajectories.
             The principal component scores can be 
             estimated through conditional expectation or via classical integration. 
             For the latter, one can choose a shrinkage method for estimated scores.
}
   
\usage{
FPCA(y,tt, p)
}

\arguments{
    \item{y}{Matrix of measurements. Each row corresponds to one subject. NAs are used in the case of sparse irrgular data 
             and regular but missing data. See help(exampleData0), help(exampleData1) for details.}
    \item{tt}{Matrix of corresponding measurement times for each measurement.} 
    \item{p}{a list obtained from setOptions(), sets the rest of arguments for FPCA. Type help(setOptions) for details.}

}
\value{
a list that contains all returned values from FPCA, where the last element contains the names of yy
  
    To see the names for yy, type names(yy). 
    To access individual values, type getVal(yy,varname). Type help(getVal) for details.
}


\references{
Yao, F., M&uuml;ller, H.G., Clifford, A.J., Dueker, S.R., Follett, J. Lin, Y., Buchholz, B. A., Vogel, J.S. (2003). Shrinkage estimation for functional component scores with application to the population kinetics of plasma folate. Biometrics 59 676-685.(pdf)

Yao, F., M&uuml;ller, H.G., Wang, J.L. (2005). Functional data analysis for sparse longitudinal data. Journal of the American Statistical Association 100 577-590. (pdf)

Yao, F., M&uuml;ller, H.G., Wang, J.L. (2005). Functional linear regression analysis for longitudinal data. Annals of Statistics 33 2873-2903. (pdf)

M&uuml;ller, H.G., Stadtm&uuml;ller, U. (2005). Generalized functional linear models. Annals of Statistics 33, 774-805. (pdf)

Chiou, J.M., M&uuml;ller, H.G. (2007). Diagnostics for functional regression via residual processes. Computational Statistics and Data Analysis 51 4849-4863. (pdf)

M&uuml;ller, H.G., Chiou, J.M., Leng, X. (2008). Inferring gene expression dynamics via functional regression analysis. BMC Bioinformatics 9:60 (pdf)

Peng, J., M&uuml;ller, H.G. (2008). Distant-based clustering of sparsely observed stochastic processes, with applications to online auctions. Annals of Applied Statistics 2 1056-1077. (pdf)

Tang, R., M&uuml;ller, H.G. (2008). Pairwise curve synchronization for functional data. Biometrika, 95 875-889. (pdf)

M&uuml;ller, H.G., Yao, F., (2008), Functional additive models. Journal of the American Statistical Association 103 1534-1544. (pdf)

Liu, B., M&uuml;ller, H.G., (2009). Estimating derivatives for samples of sparsely observed functions, with application to on-line auction dynamics. Journal of the American Statistical Association 104, 704-717. (pdf)

M&uuml;ller, H.G., Stadtm&uuml;ller, U., Yao, F., (2006). Functional variance processes. Journal of the American Statistical Association 101, 1007-1018. (pdf)

 

}


\author{Bitao Liu}

\keyword{methods}

\examples{

##load data
data(exampleData0)
y = exampleData0$y[1:50,]
tt = exampleData0$tt[1:50,]

p = setOptions(regular = 0, selection_k = "BIC1", screePlot = 1, designPlot = 1, numBins = 0, verbose = "on");
yy = FPCA(y,tt, p)
##extract some of the results for plots
out1 = getVal(yy,"out1") #vector of time points for mu, phi and y_pred
mu = getVal(yy,"mu")     #estimated mean function that is evaluated at out1
out21 = getVal(yy,"out21") #vector of time points for covariance function
xcovfit = getVal(yy, "xcovfit") #fitted covariance surface evaluated at out21 for both x and y direction
xcov = getVal(yy, "xcov")       #estimated smooth covariance evaluated at out21 for both x and y direction
phi = getVal(yy, "phi")         #estimated eigenfunctions evaluated at out1
no_opt = getVal(yy, "no_opt")   #number of estimated optimal FPC scores


#fitted covariance surface
persp(out21,out21,xcovfit,col = "lightblue",border = TRUE,
ticktype = "detailed", expand = .4, theta = 325, phi = 20, ltheta = 120,
xlab = "t", ylab = "t", main = "Fitted covariance surface", zlab = "",shade = .3)


#smooth covariance surface
persp(out21,out21,xcov,col = "lightblue",border = TRUE,ticktype = "detailed", 
expand = .4, theta = 325, phi = 20, ltheta = 120,xlab = "t", ylab = "t", 
main = "Smooth covariance surface", zlab = "",shade = .3)


#Plot observed, and predicted curves 
ypred = getVal(yy, "y_pred") 

#plot 4 randomly selected curves

n = dim(tt)[1]

rID = sample(1:n, 4, replace = FALSE)
xlim = c(min(out1),max(out1))
tmp1 = as.vector(y[!is.na(y)])
tmp2 = as.vector(ypred)
ylim = c(min(tmp1,tmp2), max(tmp1,tmp2))
par(mfrow = c(2,2))
for(i in rID){
   plot(tt[i,], y[i,], pch = 8, type = "p", xlab = "t", ylab = "X(t)", 
   main = paste("Subject", i), col= "blue", xlim = xlim, ylim = ylim)
   lines(out1, ypred[i,], col = "red")
   legend("bottomright", lty = c(1,NA), pch = c(NA,8), legend = c("fitted", "obs"), inset = .01, 
   col = c("red", "blue"),box.col = "white",cex = .9)
}


}
