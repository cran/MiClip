\name{MiClip.enriched}
\alias{MiClip.enriched}
\title{ Identify enriched bins }
\description{
  This function implements the firstHMM and tries to identify enriched bins within CLIP clusters.
}
\usage{MiClip.enriched(mic,quiet=FALSE)
}
\arguments{
  \item{mic}{
    \code{mic} is an ojbect of class "MiClip" returned by \code{MiClip.read}.
  }
  \item{quiet}{
    Whether the intermediate messages should be printed.
  }
}
\details{
The function \code{MiClip.enriched} will first divide each cluster into bins of length of \code{step} bp and then calculate the average tag coverage in each bin. Then it employs HMM algorithm and Viterbi algorithm to infer enriched bins. The output is stored in \code{enriched}.
}
\value{
An object of class \code{MiClip} is returned.
  \item{raw}{The raw data matrix including chromosomes, strands, positions, total read counts and mutant read counts}
  \item{max.hmm}{The maximum number of reads in a bin or on a base.}
  \item{model.cut}{The cutoff for fitting the mixture model in the second HMM.}
  \item{max.iterats}{The maximum number of iterations allowed for both HMM iterations.}
  \item{conver.cut}{The cutoff for reaching convergence}
  \item{enriched}{
    The output of the first HMM as a data frame. \code{region_id} is the id number generated for each cluster. \code{chr}, \code{strand}, \code{start} and \code{end} specify the genomic location of each bin. \code{tag} is the rounded average tag count in each bin. \code{enriched} and \code{probability} are the inference results.
  }
}
\seealso{
\code{\link{MiClip.read}}, \code{\link{MiClip.enriched}}, \code{\link{MiClip.binding}}, \code{\link{MiClip.sum}}
}
\examples{
data(Chi,package="MiClip")
test3=MiClip.enriched(test2,quiet=FALSE) # identify enriched regions
}

