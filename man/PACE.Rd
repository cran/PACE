\name{PACE}
\alias{PACE}
\title{The PACE package: summary information}
\description{This is the main program to perform Functional Principal Component 
             Analysis (FPCA) via PACE. The principal component scores can be 
             estimated through conditional expectation or via classical integration. 
             For the latter, one can choose a shrinkage method for estimated scores.
 }
\details{The program can be used for sparsely or densely sampled random trajectories and time courses,
via the Principal Analysis by Conditional Estimation (PACE) algorithm. 
PACE does not use pre-smoothing of trajectories which can be problematic if functional data are sparsely sampled 
or measurements are corrupted with noise. 
The main function is 'FPCA'. Three simulated data sets can be called by 'data(exampleData0)','data(exampleData1)' and 'data(exampleData2)', 
for sparse irregular, regular but missing, and completely balanced cases respectively.
Type 'help(exampleData0)' to see details.
Package 'akima' are used by this package.
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
\keyword{package}
    