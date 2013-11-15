\name{MiClip.read}
\alias{MiClip.read}
\title{ Read raw sequencing data }
\description{
  Read the sequencing data and form CLIP clusters by overlapping.
}
\usage{MiClip.read(mic)
}
\arguments{
  \item{mic}{
    \code{mic} is an ojbect of class "MiClip" returned by \code{MiClip}
  }
}
\details{
The function \code{MiClip.read} calls embeded perl scripts to read SAM format file and extract mutation information. Then CLIP clusters are formed from reads that can overlap by at least 1 bp. Reads that cannot be overlapped with any other reads are discarded.
}
\value{
An object of class \code{MiClip} is returned.
  \item{raw}{The raw data matrix including chromosomes, strands, positions, total read counts and mutant read counts}
  \item{max.hmm}{The maximum number of reads in a bin or on a base.}
  \item{empirical}{A parameter used in model fitting in the first HMM}
  \item{model.cut}{The cutoff for fitting the mixture model in the second HMM.}
  \item{max.iterats}{The maximum number of iterations allowed for both HMM iterations.}
  \item{conver.cut}{The cutoff for reaching convergence}
  \item{background}{The background data frame if available}
}
\seealso{
\code{\link{MiClip.read}}, \code{\link{MiClip.enriched}}, \code{\link{MiClip.binding}}, \code{\link{MiClip.sum}}
}
\examples{
test1=MiClip(file=system.file("extdata/test.sam",package="MiClip"),mut.type="Del")
test2=MiClip.read(test1) # read raw data
}

