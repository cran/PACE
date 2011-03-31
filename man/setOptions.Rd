\name{setOptions}
\alias{setOptions}
\title{Set all arguments for FPCA}
\description{Use this function to set arguments for FPCA.}
\usage{
setOptions(bwmu = 0, bwmu_gcv = 1, bwxcov = c(0,0), bwxcov_gcv = 1, 
                      ntest1 = 30, ngrid1 = 30, selection_k = "BIC1", FVE_threshold = 0.85,
                      maxk = 20, control = "auto", regular = 0, error = 1, ngrid = 51,
                      method = "CE", shrink = 0, newdata = NULL, kernel = "gauss", 
                      numBins = NULL, yname = NULL, screePlot = 0, designPlot = 0, rho = "cv",
                      verbose = "on")
}
\arguments{
Here is a list of all the arguments you can set by setOptions().
 \item{bwmu}{scalar bwmu>=0, bandwidth for mean curve mu,
                    0:  use CV or GCV      [Default]
                     bwmu(>0): user-specified bandwidth.}
 \item{bwmu_gcv}{Only for choice bwmu = 0,
                   0: CV   (may be time-consuming)
                   1: GCV  (faster)                            [Default]}
 \item{bwxcov}{ bandwidths for covariance surface cov(X(t),X(s))
                   bwxcov(1)==0 & bwxcov(2)==0: use CV or GCV for automatic
                   selection.                                  [Default]
                   bwxcov(1)>0 & bwxcov(2)>0: user-specified bandwidths.}
 \item{bwxcov_gcv}{If setting bwxcov = [0 0], automatic bandwidth selection
                   by CV or GCV choices
                   0: CV method.
                   1: GCV method (faster)                      [Default]}
 \item{ ntest1}{integer(<=n), number of curves used for CV when choosing
                   bandwidths for smoothing the covariance surface. [Default is 30.]}
 \item{ ngrid1}{integer, number of support points for the covariance surface 
                   in the CV procedure.    [Default is 30.]}
 \item{  selection_k}{ the method for choosing the number of principal components K.
                    'AIC1': use AIC criterion with pseudo-likelihood of
                            measurements (marginal likelihood).
                    'AIC2': use AIC criterion with likelihood of measurements
                            conditional on estimated random coeffecients.
                    'BIC1': use BIC criterion with pseudo-likelihood of
                            measurements (marginal likelihood).
                                                                [Default]
                    'BIC2': use BIC criterion with likelihood of measurements
                            conditional on estimated random coeffecients.
                    'FVE' (fraction of variance explained) : use scree plot
                            approach to select number of principal components),
                     see "FVE_threshold" below.
                     'AIC_R': use regression AIC criterion (select number
                              of principal components based on the linear
                              relationship between predictor and response)
                     positive integer K: user-specified number of principal 
                                         components}
 \item{ FVE_threshold}{  a positive number that is between 0 and 1 [Default is 0.85.]
                       It is used with the option selection_k = 'FVE' to select
                       the number of principal components that explain at least
                       "FVE_threshold" of total variation (the fraction
                       of variance explained).}
 \item{  maxk}{ integer, the maximum number of principal components to consider
                  if using automatic methods to choose K.      [Default is 20.]}
 \item{ control}{ 'auto', Select K by minimizing AIC or BIC, or find the
                          first K such that the FVE_threshold is exceeded. [Default]
                  'look', a scree plot (FVE% Vs No. of PC) will be generated based
                          on K <= 15. User will be prompted to enter user-specified
                          K after viewing scree plot. This can be combined
                          with any setting of selection_k.}
 \item{regular}{ 0, sparse (or irregular) functional data.      [Default]
                  1, regular data with missing values
                  2, completely balanced (regular) data.}
 \item{ error}{0, no additional measurement error assumed.
                  1, additional measurement error is assumed.    [Default]}
 \item{ ngrid}{integer, number of support points in each direction of 
                  covariance surface when performing principal component 
                  analysis ( ngrid > K).                    [Default is 51.]}
 \item{method}{used for computing random effects xi_{ik}
                  'CE': conditional expectation method           [Default]
                  'IN': classical integration method
                  Note: 'CE' can be applied for sparse data or regular data, but
                        'IN' can be applied only in the case of regular data.}
 \item{shrink}{    indicator of whether applying shrinkage to estimates of random
                  coefficients (for regular data only)
                  0:  no shrinkage when method = 'CE' or error = 0 [Default]
                  1:  shrinkage when method = 'IN' and error =1, otherwise, this
                      will be re-set to 0.}
 \item{newdata}{a row vector of user-defined output time grids for
                  all curves. This corresponds to "out1" in the output argument
                  If newdata = [], then "out1" corresponds to the set of distinct
                  time points from the pooled data.
                  "newdata" is supposed to be a vector in ascending order on
                   the domain of the functions.                    [Default is []]}
 
 \item{ kernel}{a character string to define the kernel to be used in the
                  1-D or 2-D smoothing
                  kernel = 'epan'  ==> Epanechnikov kernel [Default for regular=2]
                           'rect'  ==> Rectangular kernel
                           'gauss'  ==> Gaussian kernel    [Default for regular=0 or regular=1]}
 \item{numBins}{ 0: no binning
                  a positive interger (>= 10): prebin the data with user-defined
                  number of bins. When numBins < 10, no binning will be performed.
                  []:  prebin the data with the following rule    [Default]
 
                  i) When the input data is regular = 1 or 2
                     m = max of n_i, where n_i is the number of repeated measurements
                     for ith subject.
                  ii) regular = 0
                     m = median of n_i
 
                  When m <= 20 subjects, no binning.
                  When n <= 5000 subjects and m <= 400, no binning.
                  When n <= 5000 subjects and m > 400, numBins = 400.
                  When n > 5000 subjects, compute
 
                  m* = max(20, (5000-n)*19/2250+400)
 
                  if m > m*, numBins = m*
                  if m <= m*, no binning}
 \item{yname}{a character string which denotes the name of the current
                  function to be estimated.               [Default is []]
                  .}
 \item{screePlot}{indicator of whether to create the scree plot
                  1  a scree plot will be created         
                  0  no scree plot will be created        [Default]}
 \item{ designPlot}{ indicator of whether to create the design plot
                  1  a design plot will be created
                  0  no design plot will be created       [Default]
                  Interpretation of design plot: All subareas of the
                  domain x domain support square of the covariance surface need
                  to be covered more or less homogeneously with pairs of design points.}
 \item{rho}{truncation threshold for the iterative residual that is used
                in the estimation of FPC scores. (see FPCscore.pdf under Help/ for more details)
                -1:  compute unadjusted FPC scores (as in previous PACE versions)
                >0:  user-defined choice of rho
                0:   do not set truncation threshold, but use iterative residuals for sigmanew (see in output description below)
                'cv-random':  a character string which specifies to use randomized leave-one-measurement-out CV approach to find
                       the optimal value of rho.        
                Note that this choice contains a random element and therefore the analysis is not precisely
                 replicable when running  the program twice.
                'cv': use non-randomized leave-one-measurement-out CV approach to find the optimal value of rho. [default]}
 \item{verbose}{a character string
                  'on': display diagnostic messages       [Default]
                  'off': suppress diagnostic messages }
} 

\value{
    Returns a list of arguments can be used in FPCA. See also FPCA.
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

\keyword{options}


