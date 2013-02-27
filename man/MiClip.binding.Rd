\name{MiClip.binding}
\alias{MiClip.binding}
\title{ Identify binding sites }
\description{
  This function implements the second HMM and tries to identify binding sites within enriched bins.
}
\usage{MiClip.binding(mic,quiet=FALSE)
}
\arguments{
  \item{mic}{
    \code{mic} is an ojbect of class "MiClip" returned by \code{MiClip.enriched}.
  }
  \item{quiet}{
    Whether the intermediate messages should be printed.
  }
}
\details{
The function \code{MiClip.binding} will first expand all adjacent enriched bins into single base pairs and then concatenate neighboring sites. So one cluster may contain multiple enriched segments, although this is rare. Then \code{MiClip.binding} employs HMM algorithm and Viterbi algorithm to infer true binding sites. The output is stored in \code{sites}.
}
\value{
An object of class \code{MiClip} is returned.
  \item{enriched}{
    The output of the first HMM as a data frame. \code{region_id} is the id number generated for each cluster. \code{chr}, \code{strand}, \code{start} and \code{end} specify the genomic location of each bin. \code{tag} is the rounded average tag count in each bin. \code{enriched} and \code{probability} are the inference results.
  }
  \item{sites}{
    The output of the second HMM as a data frame. \code{region_id} is the id number generated for the cluster where each base resides in. \code{sub_region_id} is the id number of the concatenated segment within enriched clusters. \code{chr}, \code{strand} and \code{pos} specify the genomic location of each base. \code{tag} is the read count on each base and \code{mutant} is the mutant count on each base. \code{sites} and \code{probability} are the inference results.
  }
  \item{clusters}{
    The summary of results for all CLIP clusters. \code{clusters} contains information of chromosome, strand, start position, end position, whether or not contains enriched bins and whether or not contains binding sites.
  }
}
\seealso{
\code{\link{MiClip.read}}, \code{\link{MiClip.enriched}}, \code{\link{MiClip.binding}}, \code{\link{MiClip.sum}}
}
\examples{
data(Chi,package="MiClip")
test4=MiClip.binding(test3,quiet=TRUE) # identify binding sites
}

