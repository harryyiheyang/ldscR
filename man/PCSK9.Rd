\name{PCSK9}
\alias{PCSK9}
\title{PCSK9}
\description{
  A list with 12 eQTL summary data in the PCK9 locsc
}
\usage{
  data(PCSK9)
}
\format{
  A list with 12 data.frame:
    \describe{
      \item{SNP}{SNP identifier.}
      \item{Zscore}{Zscore of eQTL data.}
      \item{N}{Sample size of eQTL data.}
      \item{A1}{Effect allele.}
      \item{A2}{Alternative.}
      \item{rownames(.)}{Gene ensembl or splicing identifier and tissue.}
    }
}
\source{
  eQTL summary data are provided by GTEx consortium (https://www.gtexportal.org/home/).
}
\keyword{datasets}
