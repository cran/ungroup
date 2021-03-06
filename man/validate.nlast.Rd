% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate.nlast}
\alias{validate.nlast}
\title{Check if \code{nlast} needs to be adjusted in order to accommodate 
\code{out.step}}
\usage{
validate.nlast(x, nlast, out.step)
}
\arguments{
\item{x}{Vector containing the starting values of the input intervals/bins.
For example: if we have 3 bins \code{[0,5), [5,10) and [10, 15)},
\code{x} will be defined by the vector: \code{c(0, 5, 10)}.}

\item{nlast}{Length of the last interval. In the example above \code{nlast} 
would be 5.}

\item{out.step}{Length of estimated intervals in output. 
Values between 0.1 and 1 are accepted. Default: 1.}
}
\description{
Check if \code{nlast} needs to be adjusted in order to accommodate 
\code{out.step}
}
\keyword{internal}
