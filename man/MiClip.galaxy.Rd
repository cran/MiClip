\name{MiClip.galaxy}
\alias{MiClip.galaxy}
\title{ Wrapper function for running MiClip on Galaxy }
\description{
  A wrapper function for running MiClip on Galaxy mirror. This function is not to be used by outside users.
}
\usage{MiClip.galaxy(file="",control=NULL,mut.type="T2C",step=5,max.hmm=100,
paired=F,suffix=NULL,empirical="auto",model.cut=0.2,max.iterats=20,
conver.cut=0.01)
}
\arguments{
  \item{file}{
    The file name (may include path name) of the mapped tag file. \code{file} can be only in SAM format and basespace. The package can work on both single-end and paired-end datsets.
  }
  \item{control}{
    The file name of the control experiment (if available). 
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
    In the first HMM, all clusters will be divided into bins of the same length of \code{step} bp and HMM will work to distinguish enriched bins from non-enriched ones. 
  }
  \item{max.hmm}{
    The maximum number of reads in a bin or on a base. This is used to keep calculation within the dynamic range of R. If this number is too large, probability values which are very small will become zero.
  }
  \item{empirical}{
    A parameter used in model fitting in the first HMM. Default is "auto" which lets the algorithm decides its value. It can be set to the estimated minimal number of overlapping tags for a reliable CLIP cluster  if default does not work. A higher value will lead to more conservative estimation.
  }
  \item{model.cut}{
    The cutoff for fitting the mixture model in the second HMM. It can be set to the estimated minimal proportion of mutation tags vs. total tags for a binding site to be reliable. Larger values will lead to more conservative predictions. It should be between 0 and 1.
  }
  \item{max.iterats}{
    The maximum number of iterations allowed for both HMM iterations.
  }
  \item{conver.cut}{
    The cutoff for reaching convergence
  }
}