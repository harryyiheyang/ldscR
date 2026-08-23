.new_ldsc_bootstrap_result <- function(result, bootstrap) {
  attr(result, "bootstrap") <- bootstrap
  class(result) <- c("ldsc_bootstrap_result", class(result))
  result
}

#' @export
`$.ldsc_bootstrap_result` <- function(x, name) {
  if (identical(name, "bootstrap")) {
    return(attr(x, "bootstrap", exact = TRUE))
  }

  NextMethod("$")
}
