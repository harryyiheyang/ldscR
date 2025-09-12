\name{Hapmap3_AFRLDSC}
\docType{data}
\alias{Hapmap3_AFRLDSC}
\title{African Population LD Score Reference Panel}
\description{
  LD Score reference panel for African (AFR) population, constructed using
  UK Biobank data. SNPs are ordered by chromosome and base pair position.
}
\usage{data("Hapmap3_AFRLDSC")}
\format{
  A data frame with 1,234,911 observations and 2 variables:
    \describe{
      \item{\code{SNP}}{SNP identifier (rsID)}
      \item{\code{LDSC}}{LD Score coefficient}
    }
}
\source{
  Constructed from UK Biobank AFR samples via PRScsx: \url{https://github.com/getian107/PRScsx}
}
\examples{
  data(Hapmap3_AFRLDSC)
  head(Hapmap3_AFRLDSC)
  dim(Hapmap3_AFRLDSC)
}
\keyword{datasets}
