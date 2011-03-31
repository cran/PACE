\name{getVal}
\alias{getVal}
\title{Extract results form FPCA}
\description{Use this function to extract results from FPCA.}
\usage{getVal(X, varname)}

\arguments{
\item{X}{the object returned from function FPCA. Type help(FPCA) for details.}
\item{varname}{The name of variable which needs to be extracted. Use names(yy) to see all variable names.}
}

\details{
varname can be one of the following:
\itemize{
\item[no_opt] the number of selected components.
\item[sigma] estimate of measurement error variance if error=1, while it is [] if error=0.
\item[lambda] estimated eigenvalues.
\item[phi] estimated principal component functions valued at distinct input time points with ascending order of all subjects, corresponding to out1.
\item[eigens] estimated principal component functions, valued at out21, ngrid of the pooled distinct time points with ascending order of all subjects, phi is an interpolated version of eigen at out1.
\item[xi_est] predictions of random coeffecients (PC scores) for n subjects.
\item[xi_var] Var(PC score)-Var(estimated PC score). The omega matrix in equation (7) of the paper, which is used to construct the point-wise C.I. for X_i(t).
\item[mu] estimated mean functions valued at distinct input time points (newdata = []), in ascending order from all subjects, corresponding to out1; when newdata is defined, corresponds to the time points from newdata, same as out1.
\item[bw_mu] selected bandwidth for smoothing mean curve.
\item[xcov] smoothed covariance surface (diagonal removed), corresponding to out21.
\item[bw_xcov] bandwidths for smoothing covariance surface.
\item[xcovfit] fitted covariance surface, corresponding to out21.
\item[AIC] AIC values, If AIC method is not applied, it is [].
\item[BIC] BIC values, If BIC method is not applied, it is [].
\item[FVE] fraction of variance explained.
\item[y_pred] the vector of predictions for the ith subject evaluated at time points from the output grid vector "out1".
\item[y_predOrig] y_predOrig{i} is the vector of predictions for the ith subject at the same time points as the input.
\item[out1] distinct input time points with ascending order from all subjects if newdata = []; otherwise, it is the same as newdata.
\item[out21] a grid of time points for which the smoothed covariance surface assumes values, i.e., ngrids from out1.
\item[y] if no binning is performed, same as the input y if binning is performed, a vector of measurements after binning for subject i, i = 1,...,n
\item[tt] if no binning is performed, same as the input t if binning is performed,  a vector of time points after binning for subject i, i = 1,...,n
\item[regular] if no binning is performed or regular = 2, the same as the input if binning is performed and regular = 0, it will be reset to regular = 1. In other words, after binning, the sparse and irregular case is analyzed as regular data with missings.
\item[rho_opt] if rho is set as 'cv', then rho_opt is the optimal rho obtained from the CV method, otherwise, it is the same as the input rho.
\item[sigmanew] if rho is set as >0, 0 or 'cv', then sigmanew is the iterative residual sum of squares. if rho is set as -1, then sigmanew is the same as output sigma.
}}

\value{
    Returns the values of varname.
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

\keyword{getVal}


