\name{Hapmap3_EASLDSC}
\docType{data}
\alias{Hapmap3_EASLDSC}
\title{East Asian Population LD Score Reference Panel}
\description{
  LD Score reference panel for East Asian (EAS) population, constructed using
  UK Biobank data. SNPs are ordered by chromosome and base pair position.
}
\usage{data("Hapmap3_EASLDSC")}
\format{
  A data frame with 1,031,412 observations and 2 variables:
    \describe{
      \item{\code{SNP}}{SNP identifier (rsID)}
      \item{\code{LDSC}}{LD Score coefficient}
    }
}
\source{
  Constructed from UK Biobank EAS samples via PRScsx: \url{https://github.com/getian107/PRScsx}
}
\examples{
  data(Hapmap3_EASLDSC)
  head(Hapmap3_EASLDSC)
  dim(Hapmap3_EASLDSC)
}
\keyword{datasets}

