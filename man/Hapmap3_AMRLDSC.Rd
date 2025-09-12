\name{Hapmap3_AMRLDSC}
\docType{data}
\alias{Hapmap3_AMRLDSC}
\title{Admixed American Population LD Score Reference Panel}
\description{
  LD Score reference panel for Admixed American (AMR) population, constructed
  using UK Biobank data. SNPs are ordered by chromosome and base pair position.
}
\usage{data("Hapmap3_AMRLDSC")}
\format{
  A data frame with 1,183,356 observations and 2 variables:
    \describe{
      \item{\code{SNP}}{SNP identifier (rsID)}
      \item{\code{LDSC}}{LD Score coefficient}
    }
}
\source{
  Constructed from UK Biobank AMR samples via PRScsx: \url{https://github.com/getian107/PRScsx}
}
\examples{
  data(Hapmap3_AMRLDSC)
  head(Hapmap3_AMRLDSC)
  dim(Hapmap3_AMRLDSC)
}
\keyword{datasets}
