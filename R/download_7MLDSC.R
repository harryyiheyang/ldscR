#' Download 7M SNP LD Score Reference Panels
#'
#' @description
#' Downloads LD Score reference panels computed from UK Biobank imputed data
#' containing approximately 7 million SNPs. These panels are derived from the
#' SBayesRC project which provides high-density LD reference data for multiple
#' ancestries and cross-population combinations.
#'
#' The LD reference panels include:
#' \itemize{
#'   \item Single ancestry panels: EUR, EAS, AFR computed from UK Biobank
#'         samples of respective ancestries using ~7M imputed SNPs
#'   \item Cross-population panels: EURxEAS, EURxAFR, EASxAFR computed by
#'         combining LD information across ancestry groups
#'   \item SNP information file (7m_snp.parquet) containing SNP identifiers
#'         and allele information for harmonization
#' }
#'
#' These high-density panels improve polygenic prediction accuracy compared to
#' traditional HapMap3-based references by leveraging the full spectrum of
#' common variation captured by modern imputation.
#'
#' @param ancestry Character string specifying which LD reference to download.
#'   Options: "EUR", "EAS", "AFR", "EURxEAS", "EURxAFR", "EASxAFR"
#'
#' @return Invisible NULL. The function loads data directly into the global environment:
#'   \item{SevenMillon_[ancestry]LDSC}{Data frame with LD scores for the requested ancestry}
#'   \item{SevenMillon_SNP}{Data frame with SNP information (loaded once per session)}
#'
#' @details
#' The function downloads parquet files from a remote server and loads them
#' directly into the global environment using assign(). The SNP information
#' file is only downloaded once per session to avoid redundant transfers.
#'
#' File sizes are substantial (~500MB-1GB each) due to the high SNP density.
#' Ensure adequate internet connection and disk space.
#'
#' @source
#' LD reference panels computed using SBayesRC methodology from:
#' \url{https://github.com/zhilizheng/SBayesRC}
#'
#' Based on UK Biobank imputed genotype data with ~7 million common SNPs,
#' providing enhanced LD estimation compared to sparse reference panels.
#'
#' @references
#' Zheng, Z., et al. (2024). Leveraging functional genomic annotations and
#' genome coverage to improve polygenic prediction of complex traits within
#' and between ancestries. Nature Genetics.
#'
#' @examples
#' \dontrun{
#' # Download European LD reference
#' download_7MLDSC("EUR")
#' head(SevenMillon_EURLDSC)  # LD scores now available in workspace
#'
#' # Download cross-population EUR x EAS reference
#' download_7MLDSC("EURxEAS")
#' dim(SevenMillon_EURxEASLDSC)  # Cross-population LD scores
#'
#' # Check SNP information (loaded automatically)
#' head(SevenMillon_SNP)
#'
#' # Check what's in your workspace
#' ls(pattern = "^SevenMillon_")
#' }
#'
#' @export
download_7MLDSC <- function(ancestry) {

  # Check if arrow package is available
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required but not installed. Please install it with: install.packages('arrow')")
  }

  # Validate ancestry input
  valid_ancestries <- c("EUR", "EAS", "AFR", "EURxEAS", "EURxAFR", "EASxAFR")
  if (!ancestry %in% valid_ancestries) {
    stop("ancestry must be one of: ", paste(valid_ancestries, collapse = ", "))
  }

  # Base URL for downloads
  base_url <- "https://ldsc7m.netlify.app"

  # Create temporary directory
  temp_dir <- tempfile()
  dir.create(temp_dir)

  # Ensure cleanup on exit
  on.exit(unlink(temp_dir, recursive = TRUE))

  cat("Downloading 7M LD Score reference for", ancestry, "ancestry...\n")

  # Construct file names
  ldsc_file <- paste0("7m_", tolower(gsub("x", "x", ancestry)), "ldsc.parquet")
  snp_file <- "7m_snp.parquet"

  # Download LDSC file
  ldsc_url <- file.path(base_url, ldsc_file)
  ldsc_path <- file.path(temp_dir, ldsc_file)

  cat("  Downloading", ldsc_file, "...\n")
  tryCatch({
    download.file(ldsc_url, ldsc_path, mode = "wb", quiet = TRUE)
  }, error = function(e) {
    stop("Failed to download ", ldsc_file, ": ", e$message)
  })

  # Check if SNP file already exists in working directory
  snp_exists <- file.exists(snp_file)

  # Check if SevenMillon_SNP already exists in global environment
  sevenm_snp_exists <- exists("SevenMillon_SNP", envir = .GlobalEnv)

  if (!snp_exists && !sevenm_snp_exists) {
    # Download SNP file
    snp_url <- file.path(base_url, snp_file)
    snp_path <- file.path(temp_dir, snp_file)

    cat("  Downloading", snp_file, "...\n")
    tryCatch({
      download.file(snp_url, snp_path, mode = "wb", quiet = TRUE)
    }, error = function(e) {
      stop("Failed to download ", snp_file, ": ", e$message)
    })

    # Load SNP data and assign to global environment
    cat("  Loading SNP information...\n")
    snp_data <- arrow::read_parquet(snp_path)
    assign("SevenMillon_SNP", snp_data, envir = .GlobalEnv)

  } else if (snp_exists && !sevenm_snp_exists) {
    cat("  SNP file already exists, loading from working directory...\n")
    snp_data <- arrow::read_parquet(snp_file)
    assign("SevenMillon_SNP", snp_data, envir = .GlobalEnv)

  } else {
    cat("  SevenMillon_SNP already loaded in workspace\n")
  }

  # Load LDSC data and assign to global environment
  cat("  Loading LD scores...\n")
  ldsc_data <- arrow::read_parquet(ldsc_path)

  # Create variable name
  var_name <- paste0("SevenMillon_", ancestry, "LDSC")
  assign(var_name, ldsc_data, envir = .GlobalEnv)

  # Report dimensions and what was loaded
  cat("  Loaded", nrow(ldsc_data), "SNPs with LD scores as", var_name, "\n")
  if (exists("SevenMillon_SNP", envir = .GlobalEnv)) {
    cat("  SNP information available as SevenMillon_SNP with", nrow(get("SevenMillon_SNP", envir = .GlobalEnv)), "SNPs\n")
  }

  cat("Download complete! Data loaded to workspace.\n")

  # Return invisibly
  invisible(NULL)
}
