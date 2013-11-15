\name{MiClip}
\alias{MiClip}
\title{ A Model-based Approach to Identify Binding Sites in CLIP-Seq Data }
\description{
  Construct a "MiClip" class object for following analysis.
}
\usage{MiClip(file="",mut.type="T2C",step=5,max.hmm=100,paired=F,suffix=NULL,
  empirical="auto",model.cut=0.2,max.iterats=20,conver.cut=0.01,background=NULL)
}
\arguments{
  \item{file}{
    The file name (may include path name) of the mapped tag file. \code{file} must be in SAM format and basespace but either single-end or paired-end mode.
  }
  \item{mut.type}{
    The marker mutation for the CLIP-Seq experiment, separated by ",", e.g. "T2C", "T2C,T2A" or "T2C,Ins,Del". "T2C" denotes T-to-C substitution, "Ins" denotes insertion of any length and "Del" denotes deletion of any length. The default is "T2C". If \code{mut.type} is set to "all", all kinds of mutations are included as marker mutation.
  }
  \item{paired}{
    Whether the sequencing data is paired-end. Default is FALSE. 
  }
  \item{suffix}{
    The suffix of the paired-end read data. This is a vector which contains the suffix of the names of forward reads and backward reads. For example, if the mate pairs in the SAM file are named as "1_2_100708_26_788_F3", "1_2_100708_26_788_F5-RNA", etc, \code{suffix} can either be \code{c("F3","F5-RNA")} or \code{c("_F3","_F5-RNA")}. Default is NULL and will be set automatically to \code{c("1","2")} if \code{paired} is TRUE but \code{suffix} is not set.
  }
  \item{step}{
    In the first HMM, all clusters will be divided into bins of the same length of \code{step} bp and HMM will work to distinguish enriched bins from non-enriched ones. Default is 5 and for larger dataset(>20M mapped reads) it is better to set \code{step} to a value between 10-15.
  }
  \item{empirical}{
    Used to help model fitting in the first HMM. Default is "auto" which lets the algorithm decides its value. It can be set to the estimated minimal number of overlapping tags for a reliable enriched CLIP cluster  if default does not work. A higher value will lead to more conservative estimation.
  }
  \item{model.cut}{
    The cutoff for fitting the mixture model in the second HMM. It can be set to the estimated minimal proportion of mutation tags vs. total tags for a binding site to be reliable. Larger values will lead to more conservative predictions. It should be between 0 and 1 and the default is 0.2.
  }
  \item{max.hmm}{
    The maximum number of tag counts in a bin or on a base. This is used to keep calculation within the dynamic range of R. If this number is too large, probability values which are very small will become zero. Default is 100.
  }
  \item{max.iterats}{
    The maximum number of iterations allowed for both HMM iterations. Default is 20.
  }
  \item{conver.cut}{
    The cutoff for reaching convergence. Default is 0.01.
  }
  \item{background}{
    An optional data frame with preprocessed count data in the control experiment for normalization purpose. It must be of the same format as the \code{raw} data frame in the output of \code{MiClip.read}. The first 4 columns must be the same as \code{raw} and the corresponding columns for total tag count should be replaced with tag count in the control condition. We suggest running \code{MiClip} without the background experiment first, constructing a data frame according to the \code{raw} data frame, and then running \code{MiClip} again with that data frame.
  }
}
\details{
The function \code{MiClip} takes all the necessary parameters for calculation and constructs the initial \code{MiClip} class object.
}
\value{
An object of class \code{MiClip} is returned.
  \item{file}{The file name (may include path name) of the alignment file.}
  \item{mut.type}{The type of mutation wanted.}
  \item{paired}{Whether the sequencing data is paired-end.}
  \item{suffix}{The suffix of the paired-end read data.}
  \item{step}{Bin length.}
  \item{max.hmm}{The maximum number of reads in a bin or on a base.}
  \item{empirical}{A parameter used in model fitting in the first HMM}
  \item{model.cut}{The cutoff for fitting the mixture model in the second HMM.}
  \item{max.iterats}{The maximum number of iterations allowed for both HMM iterations.}
  \item{conver.cut}{The cutoff for reaching convergence}
  \item{background}{background data frame}
}
\seealso{
\code{\link{MiClip.read}}, \code{\link{MiClip.enriched}}, \code{\link{MiClip.binding}}, \code{\link{MiClip.sum}}
}
\examples{
test1=MiClip(file=system.file("extdata/test.sam",package="MiClip"),mut.type="Del")
# for paired-end data
# test=MiClip(file="test.sam",paired=TRUE,suffix=c("F3","F5-RNA"))
}

