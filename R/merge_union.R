#' Filter and Align GWAS Data to a Reference Panel
#'
#' The \code{merge_union} function processes a list of GWAS summary statistics data frames, harmonizes alleles according to a reference panel, removes duplicates, and aligns data to the union of SNPs. It's used to prepare data for further analysis such as LDSC.
#'
#' @param gwas_data_list A list of data.frames where each data.frame contains GWAS summary statistics for a trait. Each data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param ref_panel A data.frame containing the reference panel data. It must include columns for SNP, A1, and A2.
#'
#' @return A list of data.frames, each corresponding to an input GWAS summary statistics data frame, but filtered, harmonized, and aligned to the union of SNPs SNPs found across all data frames.
#'
#' @examples
#' # Assuming GWAS_List and ref_panel are already defined:
#' GWAS_List <- merge_union(GWAS_List, ref_panel)
#'
#' @details The function performs several key steps: adjusting alleles according to a reference panel, removing duplicate SNPs, and aligning all GWAS data frames to the union of SNPs SNPs. This is often a necessary preprocessing step before performing genetic correlation and heritability analyses.
#' @importFrom data.table setDT data.table
#' @export

merge_union <- function(gwas_data_list, ref_panel, missing.thres = 0.8, allele_match=T) {
  print("Adjusting effect allele according to reference panel...")
  p <- length(gwas_data_list)
  NAM <- names(gwas_data_list)
  if(allele_match==T){
  for (i in 1:p) {
    A <- gwas_data_list[[i]]
    A <- setDT(A)
    A <- allele_harmonise(ref_panel = ref_panel, gwas_data = A)
    A <- A[, .(SNP, Zscore, N)]
    gwas_data_list[[i]] <- A
  }
  }
  print("Finding union of SNPs...")
  snp_sets <- lapply(gwas_data_list, function(df) {
    return(as.character(df$SNP))
  })

  union_snps <- Reduce(union, snp_sets)

  print("Aligning data to union of SNPs and ordering...")
  gwas_data_union_aligned <- lapply(gwas_data_list, function(df, union_snps) {
    # Create a template data frame with all SNPs and merge it with the actual data
    template_df <- data.table(SNP = union_snps)
    df_union <- merge(template_df, df, by = "SNP", all.x = TRUE)
    df_union <- df_union[order(df_union$SNP), ]
    return(df_union)
  }, union_snps)

  print("Removing SNPs and Exposures with large number of missing values")
  gwas_data_union_aligned <- missing_value(gwas_data_list = gwas_data_union_aligned, missing.thres = missing.thres)

  print("Filtering complete.")
  return(gwas_data_union_aligned)
}

