filterss <- function(gwas_data_list, ref_panel) {
  print("Adjusting effect allele according to reference panel...")
  p <- length(gwas_data_list)
  for (i in 1:p) {
    A <- gwas_data_list[[i]]
    A <- allele_harmonise(ref_panel = ref_panel, gwas_data = A)
    gwas_data_list[[i]] <- A
  }

  gwas_data_list <- lapply(gwas_data_list, function(df) {
    df <- df[!duplicated(df$SNP), ]
    return(df)
  })

  print("Finding common SNPs...")
  snp_sets <- lapply(gwas_data_list, function(df) {
    return(as.character(df$SNP))
  })

  common_snps <- Reduce(intersect, snp_sets)

  print("Aligning data to common SNPs and ordering...")
  gwas_data_common_aligned <- lapply(gwas_data_list, function(df) {
    df_common <- df[df$SNP %in% common_snps, ]
    df_common <- df_common[order(df_common$SNP), ]
    rownames(df_common) <- 1:nrow(df_common)
    return(df_common)
  })

  print("Filtering complete.")
  return(gwas_data_common_aligned)
}
