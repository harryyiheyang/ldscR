#' Download and Load Precomputed LDSC Reference Data
#'
#' This function downloads `.rda` files (e.g., EURLDSC, AFRLDSC, EASLDSC, SNPInfo)
#' from an external Netlify server and optionally loads the data into the global environment.
#'
#' @param file A character string specifying which LDSC reference file to download.
#'             Must be one of `"eurldsc"`, `"afrldsc"`, `"easldsc"` or `"snpinfo"`.
#' @param load_into_env Logical; if \code{TRUE} (default), the downloaded file will be loaded into \code{.GlobalEnv}.
#'
#' @return Invisibly returns the path to the downloaded `.rda` file.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download and load EURLDSC into current session
#' download_ldsc_data("eurldsc")
#'
#' # Just download SNPInfo without loading
#' path <- download_ldsc_data("snpinfo", load_into_env = FALSE)
#' }
download_ldsc_data <- function(file = c("eurldsc", "afrldsc", "easldsc","snpinfo"),
                               load_into_env = TRUE) {
  file <- match.arg(file)
  base_url <- "https://ldsc7m.netlify.app/"
  url <- paste0(base_url, file, ".rda")
  dest <- file.path(tempdir(), paste0(file, ".rda"))

  message("Downloading: ", url)
  utils::download.file(url, dest, mode = "wb")

  if (load_into_env) {
    message("Loading into environment: ", file)
    load(dest, envir = .GlobalEnv)
  }

  invisible(dest)
}
