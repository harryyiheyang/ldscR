\name{Hapmap3}
\docType{data}
\alias{Hapmap3}
\title{HapMap3 SNP Information for Allele Harmonization}
\description{
  SNP information dataset containing approximately 1.7 million SNPs in common
  with either HapMap3 or the UK Biobank, used for allele harmonization in
  genetic analyses.
}
\usage{data("Hapmap3")}
\format{
  A data frame with 1,664,852 observations and 3 variables:
    \describe{
      \item{\code{SNP}}{SNP identifier (rsID)}
      \item{\code{A1}}{Effect allele}
      \item{\code{A2}}{Reference allele}
    }
}
\source{
  Derived from bigsnpr package: \url{https://privefl.github.io/bigsnpr/reference/download_1000G.html}
}
\examples{
  data(Hapmap3)
  head(Hapmap3)
  dim(Hapmap3)
}
\keyword{datasets}
