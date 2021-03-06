\name{exampleData1}
\docType{data}
\alias{exampleData1}
\title{An example for regular data with missings}
\description{A simulated dataset}
\format{
example is a list with eight components (in the given order):
\describe{
\item{ncohort}{Sample Size.}
\item{ni}{Number of measurement for each subject.}
\item{maxNi}{Maximum number of measurements.}
\item{y}{A matrix (row=ncohort, col=maxNi) of measurements. Missing values are set to be NA.}
\item{tt}{A matrix (row=ncohort, col=maxNi) of measurement times. If the measurement is missing, this time point is set to be NA.}
\item{lambda}{True eigenvalues.}
\item{xi}{True eigen-components}
\item{sigma2}{Variance of the measurement error}
 }
}

\details{This is a simulated regular data set with missings. Maximum number of measurements is 20, minimum is 1. Missing points are chosen randomly .
    true mean function: t+sin(t),
    1st eigenfunction: -sqrt(0.2)*cos(2*t*pi/10),
    2nd eigenfuncton: sqrt0.2)*sin(2*t*pi/10),
    measurement error is distributed as N(0,Sigma2).}
\usage{data(exampleData1)}
\keyword{datasets}
