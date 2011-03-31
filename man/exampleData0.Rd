\name{exampleData0}
\docType{data}
\alias{exampleData0}
\title{An example for sparse irregular data}
\description{A simulated dataset as an example}
\format{
example is a list with eight components (in the given order):
\describe{
\item{ncohort}{Sample Size.}
\item{ni}{Number of measurements for each subject.}
\item{maxNi}{Maximum number of measurements.}
\item{y}{A matrix (row=ncohort, col=maxNi) of measurements. For row i, the first Ni[i] components are the measurements, and the rest are set to be NA.}
\item{tt}{A matrix (row=ncohort, col=maxNi) of measurement times. For row i, the first Ni[i] components are the time points, and the rest are set to be NA.}
\item{lambda}{True eigenvalues.}
\item{xi}{True eigen-components}
\item{sigma2}{Variance of the measurement error}
 }
}

\details{This is a simulated sparse irregular data set. Maximum number of measurements is 4, minimum is 1. The measurement times 
are uniformly distrbuted on [0,10].
    true mean function: t+sin(t),
    1st eigenfunction: -sqrt(0.2)*cos(2*t*pi/10),
    2nd eigenfunction: sqrt0.2)*sin(2*t*pi/10),
    measurement error is distributed as N(0,Sigma2).}
\usage{data(exampleData0)}
\keyword{datasets}
