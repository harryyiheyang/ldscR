\name{Hapmap3_EURLDSC}
\docType{data}
\alias{Hapmap3_EURLDSC}
\title{European Population LD Score Reference Panel}
\description{
  LD Score reference panel for European (EUR) population, constructed using
  UK Biobank data. SNPs are ordered by chromosome and base pair position.
}
\usage{data("Hapmap3_EURLDSC")}
\format{
  A data frame with 1,117,425 observations and 2 variables:
    \describe{
      \item{\code{SNP}}{SNP identifier (rsID)}
      \item{\code{LDSC}}{LD Score coefficient}
    }
}
\source{
  Constructed from UK Biobank EUR samples via PRScsx: \url{https://github.com/getian107/PRScsx}
}
\examples{
  data(Hapmap3_EURLDSC)
  head(Hapmap3_EURLDSC)
  dim(Hapmap3_EURLDSC)
}
\keyword{datasets}
