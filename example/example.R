library(devtools)
data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist=readRDS(temp_file)
unlink(temp_file)
load_all()
data("Hapmap3")
gwaslist=filter_align(gwas_data_list=gwaslist,ref_panel=Hapmap3[,c("SNP","A1","A2")])
fit1=ldsc.univ(gwas=gwaslist$driving,LDSC=Hapmap3_EURLDSC,nblock=500,sampling.time=200)
fit2=ldsc.univ(gwas=gwaslist$computer,LDSC=Hapmap3_EURLDSC,nblock=500,sampling.time=200)
fit3=ldsc.bicov(gwas1=gwaslist$driving,gwas2=gwaslist$computer,LDSC=Hapmap3_EURLDSC,nblock=500,sampling.time=200,h21=fit1$h2,h22=fit2$h2)
fit1
fit2
fit3
gcov_matrix=matrix(c(fit1$h2,fit3$gcov,fit3$gcov,fit2$h2))
ecov_matrix=matrix(c(fit1$intercept,fit3$ecov,fit3$ecov,fit2$intercept))
