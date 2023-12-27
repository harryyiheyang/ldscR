missing_value=function(gwas_data_list, missing.thres = 0.8){
p = length(gwas_data_list)
snpvec=seq(missing.thres,0.95,by=(0.95-missing.thres)/3)
snpvec=sort(snpvec, decreasing = T)
M=matrix(0, nrow(gwas_data_list[[1]]), p)

for(i in 1:p){
M[,i] = gwas_data_list[[i]]$Zscore
}
rownames(M)=gwas_data_list[[1]]$SNP
colnames(M)=names(gwas_data_list)

for(i in 1:length(snpvec)){
M=demissing(M, rowfirst = F, rowthres = snpvec[i], colthres = snpvec[i])
}

snplist=rownames(M)
exposure_remove=setdiff(names(gwas_data_list),colnames(M))
if(length(exposure_remove)>0){
  gwas_data_list <- gwas_data_list[!(names(gwas_data_list) %in% exposure_remove)]
}

updated_gwas_data_list <- lapply(gwas_data_list, function(df) {
  df <- df[df$SNP %in% snplist, ]
  df$Zscore[is.na(df$Zscore)] <- 0
  df$N[is.na(df$N)] <- 1
  return(df)
})

return(A=list(gwas_data_list=updated_gwas_data_list,exposure_remove=exposure_remove))
}
