library(devtools)
document()
data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist=readRDS(temp_file)
unlink(temp_file)
download_ldsc_data("snpinfo")
gwaslist=filter_align(gwas_data_list=gwaslist,ref_panel=SNPInfo)
download_ldsc_data("eurldsc")
fitldsc=ldscR(GWAS_List=gwaslist,LDSC=EURLDSC,estimate_SE=T)
fitldsc$GCovEst
fitldsc$GCovSE
fitldsc$GCorEst
fitldsc$GCorSE
