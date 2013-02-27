\name{MiClip.snp}
\alias{MiClip.snp}
\title{ Distinguish possible SNPs from high confidence CLIP binding sites }
\description{
  Read the alignment data of the control experiments (e.g. RNA with no cross-linking) and mark those high confidence CLIP bindings sites which might be SNPs.
}
\usage{MiClip.snp(mic,file="",mut.type="T2C",paired=F,suffix=NULL)
}
\arguments{
  \item{mic}{
    \code{mic} is an ojbect of class "MiClip" returned by \code{MiClip.binding}
  }
  \item{file}{
    \code{file} is the alignment file of the control experiments. 
  }
  \item{mut.type}{
    \code{mut.type} must be the same as the one used for \code{MiClip} when running the crosslinked sample.
  }
  \item{paired}{
    \code{paired} must be the same as the one used for \code{MiClip} when running the crosslinked sample.
  }
  \item{suffix}{
    \code{suffix} must be the same as the one used for \code{MiClip} when running the crosslinked sample.
  }
}
\details{
The function \code{MiClip.snp} is devised to distinguish possible SNPs from inferred high-confidence crosslinking sites. Optionally, users can add additional quality control steps before alignment of the control data to the reference genome. \code{MiClip.snp} takes the alignment file of the control condition and looks for the same mutations as in the treatment sample. A null hypethesis is tested on each mutant site by \code{MiClip.snp} in order to extracts possible SNPs. Then the binding sites inferred by \code{MiClip.binding} are screened for possible SNPs (mutant sites that are not inferred as binding sites are ignored). A column will be added to \code{sites} in the final output specifying whether a binding site could actually be a SNP. And another column will be added to \code{clusters} in the final output specifying those clusters, at least one of whose binding sites could be a SNP.    
}
\value{
An object of class \code{MiClip} is returned.
  \item{enriched}{
    The same as the one returned by \code{MiClip.binding}.
  }
  \item{sites}{
    The same as the one returned by \code{MiClip.binding}.
  }
  \item{clusters}{
    The same as the one returned by \code{MiClip.binding}.
  }
  \item{snps}{
    This data frame is added by \code{MiClip.snp}. It contains information of the possible SNP sites extracted from the control experiment file.
  }
}
\seealso{
\code{\link{MiClip.read}}, \code{\link{MiClip.enriched}}, \code{\link{MiClip.binding}}, \code{\link{MiClip.sum}}
}
\examples{
data(Chi,package="MiClip")
test5=MiClip.snp(test4,file=system.file("extdata/snp.sam",package="MiClip"),mut.type="Del") # mark possible SNPs
}

