\name{MiClip.sum}
\alias{MiClip.sum}
\title{ Summary of MiClip Inference Results }
\description{
  This summary function computes simple statistics for the results produced by \code{MiClip.binding}.
}
\usage{
MiClip.sum(mic,...)
}
\arguments{
  \item{mic}{
    \code{mic} is an ojbect of class "\code{MiClip}" returned by \code{MiClip.enriched} or \code{MiClip.binding}.
  }
  \item{...}{
    further arguments passed to or from other methods.
  }
}
\details{
  This function will compute summary statistics only if \code{mic} is generated from \code{MiClip.binding}.
}
\seealso{
\code{\link{MiClip.read}}, \code{\link{MiClip.enriched}}, \code{\link{MiClip.binding}}, \code{\link{MiClip.sum}}
}
\examples{
data(Chi,package="MiClip")
MiClip.sum(test4) # print summary
}

