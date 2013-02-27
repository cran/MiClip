\name{MiClip.adaptor}
\alias{MiClip.adaptor}
\title{ Trim 3' adaptor }
\description{
  This helper function will remove 3' adaptors from raw reads in the sequence file.
}
\usage{MiClip.adaptor(file="",format="fastq",adaptor="",min=15,mismatch=0.4)
}
\arguments{
  \item{file}{
    The filename (including full path name) of the sequencing file.
  }
  \item{format}{
    The format of the sequencing file. It can be either "fastq" or "fasta". Also the raw sequencing file must be in basespace.
  }
  \item{adaptor}{
    The adaptor sequence, for example "TCGTATGCCGTCTTCTGCTTG". "N" is allowed, and it is case insensitive. So "TCGTNNGCCGTCttcnncttg" is also ok.
  }
  \item{min}{
    After trimming, if a sequence is shorter than \code{min}, it will be tossed away.
  }
  \item{mismatch}{
    The maximum proportion of mismatches allowed when aligning adaptor sequence to the 3' end. 
  }
}
\details{
  This function is a wrapper function of a perl script. It trimms a full or partial 3' adaptor from each sequencing read and generates a new file in the same folder of the original sequencing file. It can only work on single-end reads now.
  For example, if \code{adaptor} is "TCGTATGCCGTCTTCTGCTTG", \code{min} is 15 and \code{mismatch} is 0.25, "NNTGGAGGCCGGACGCTTCCNAAANNNGTATGTCGT" will be trimmed down to "NNTGGAGGCCGGACGCTTCCNAAAN". There is only one mismatch in the partial adaptor sequence "NNGTATGTCGT" and 1/11<0.25, so this part will be trimmed from the short read.
  The adaptor at the 5' end usually won't be sequenced. Even if part of the 5' end adaptor is sequenced, such cases are usually rare. So 5' end adaptor is not considered in this function. If the user would like to remove 5' end adaptor too, please refer to other specialized adaptor removing algorithm. 
}
\examples{
library("MiClip")

MiClip.adaptor(file=system.file("extdata/test.fastq",package="MiClip"),
  adaptor="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC")
}

