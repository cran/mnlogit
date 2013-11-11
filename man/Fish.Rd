\name{Fish}
\docType{data}
\alias{Fish}
\alias{Fishing}
\title{Choice of Fishing Mode }
\description{
  A data.frame containing data on choices of recreational fishing mode. 
  Data may depend on both the individual and the alternative.

 \emph{number of observations} :  1182

 \emph{country} :  United States

}
\usage{data(Fish)}
\format{A dataframe containing :
  \describe{
    \item{mode}{- recreation mode choice, one of : beach, pier, boat and charter}
    \item{price}{- price for a mode for an individual}
    \item{catch}{- fish catch rate for a mode for an individual}
    \item{income}{- monthly income of an individual}
    \item{chid}{- chooser ID: serial number of the individual}
  }
}
\source{
     Data taken from the R package \pkg{mlogit} by Yves Croissant, which lists the
     source as:
 
     Herriges, J. A.  and C. L.  Kling (1999) \dQuote{Nonlinear Income
     Effects in Random Utility Models}, \emph{Review of Economics and
     Statistics}, \bold{81}, 62-72.
}
\references{
    Cameron, A.C.  and P.K.  Trivedi (2005) \emph{Microeconometrics :
    methods and applications}, Cambridge, pp. 463--466, 486 and
    491--495.
}

\keyword{datasets}
