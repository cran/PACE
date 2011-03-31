\name{exampleData2}
\docType{data}
\alias{exampleData2}
\title{An example for completely balanced data}
\description{A simulated dataset}
\format{
example is a list with seven components (in the given order):
\describe{
\item{ncohort}{Sample Size.}
\item{maxNi}{Maximum number of measurements.}
\item{y}{A matrix (row=ncohort, col=maxNi) of measurements.}
\item{tt}{A matrix (row=ncohort, col=maxNi) of measurement times.}
\item{lambda}{True eigenvalues.}
\item{xi}{True eigen-components}
\item{sigma2}{Variance of the measurement error}
 }
}

\details{This is a simulated completely balanced data set.
    true mean function: t+sin(t),
    1st eigenfunction: -sqrt(0.2)*cos(2*t*pi/10),
    2nd eigenfuncton: sqrt0.2)*sin(2*t*pi/10),
    measurement error is distributed as N(0,Sigma2).}
\usage{data(exampleData2)}
\keyword{datasets}
