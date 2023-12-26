\name{AFRLDSC}
\alias{AFRLDSC}
\title{African Population LDSC Data}
\description{
  LDSC data for African populations. It includes LD scores calculated from the sample LD matrix and POET-estimated LD scores.
}
\usage{
  data(AFRLDSC)
}
\format{
  A data frame with 2 variables:
    \describe{
      \item{SNP}{SNP identifier.}
      \item{LDSC}{LD score calculated from the sample LD matrix.}
    }
}
\source{
  LD scores derived from UK Biobank samples using PRScsx (https://github.com/getian107/PRScsx).
}
\keyword{datasets}
