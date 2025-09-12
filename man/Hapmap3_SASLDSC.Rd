\name{Hapmap3_SASLDSC}
\docType{data}
\alias{Hapmap3_SASLDSC}
\title{South Asian Population LD Score Reference Panel}
\description{
  LD Score reference panel for South Asian (SAS) population, constructed using
  UK Biobank data. SNPs are ordered by chromosome and base pair position.
}
\usage{data("Hapmap3_SASLDSC")}
\format{
  A data frame with 1,133,674 observations and 2 variables:
    \describe{
      \item{\code{SNP}}{SNP identifier (rsID)}
      \item{\code{LDSC}}{LD Score coefficient}
    }
}
\source{
  Constructed from UK Biobank SAS samples via PRScsx: \url{https://github.com/getian107/PRScsx}
}
\examples{
  data(Hapmap3_SASLDSC)
  head(Hapmap3_SASLDSC)
  dim(Hapmap3_SASLDSC)
}
\keyword{datasets}
