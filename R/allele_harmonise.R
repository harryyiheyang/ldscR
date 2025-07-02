#' allele_harmonise: Function to harmonize allele coding between a reference panel and GWAS/eQTL summary data
#'
#' This function harmonizes the allele coding between an LD reference panel and GWAS/eQTL summary data to ensure consistency in the effect allele (A1) and the other allele (A2). The Z-scores in the GWAS/eQTL summary data are adjusted accordingly if the alleles are flipped.
#'
#' @param ref_panel A data.frame or data.table representing an LD reference panel. It must contain columns "SNP", "A1" (effect allele), and "A2" (other allele).
#' @param gwas_data A data.frame or data.table representing GWAS/eQTL summary data. It must contain columns "SNP", "A1", "A2", and "Zscore".
#'
#' @return A data.table with harmonized alleles and adjusted Z-scores. The output contains the harmonized alleles from the reference panel and the original alleles from the GWAS/eQTL data as separate columns.
#' @importFrom data.table setDT setkey setnames
#' @export
#'
allele_harmonise <- function(ref_panel, gwas_data) {

  # Make sure the reference panel has columns SNP, A1 and A2
  stopifnot(all(c("SNP", "A1", "A2") %in% colnames(ref_panel)))

  # Make sure the GWAS data has columns SNP, A1, A2, and Tstat
  stopifnot(all(c("SNP", "A1", "A2", "Zscore") %in% colnames(gwas_data)))

  # Merge the reference panel with the GWAS data
  merged_data <- merge(gwas_data, ref_panel, by = "SNP")

  # Remove any SNPs where A1 and A2 do not match between the reference panel and the GWAS data
  ind=which((toupper(merged_data$A1.x) == toupper(merged_data$A1.y) & toupper(merged_data$A2.x) == toupper(merged_data$A2.y))
            |(toupper(merged_data$A1.x) == toupper(merged_data$A2.y) & toupper(merged_data$A2.x) == toupper(merged_data$A1.y)))

  merged_data=merged_data[ind,]

  # Multiply Tstat by -1 if A1 and A2 are opposed between the reference panel and the GWAS data
  merged_data$Zscore <- ifelse(toupper(merged_data$A1.x) == toupper(merged_data$A2.y) & toupper(merged_data$A2.x) == toupper(merged_data$A1.y), -1 * merged_data$Zscore, merged_data$Zscore)

  colnames(merged_data)[which(names(merged_data) == "A1.y")] <- "A1"
  colnames(merged_data)[which(names(merged_data) == "A2.y")] <- "A2"
  colnames(merged_data)[which(names(merged_data) == "A1.x")] <- "A1.orginal"
  colnames(merged_data)[which(names(merged_data) == "A2.x")] <- "A2.orginal"

  return(merged_data)
}
