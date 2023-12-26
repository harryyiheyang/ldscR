# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(ldscR)
#test_check("ldscR")

url <- "http://tinyurl.com/5dszcsfz"
dest_file <- "GWAS_List.rds"
download.file(url, destfile = dest_file, mode = "wb")
gwas_list <- readRDS(dest_file)
file.remove(dest_file)
GWAS_AFR=gwas_list[c(1,5,9,13)]
GWAS_EAS=gwas_list[c(1,5,9,13)+1]
GWAS_EUR=gwas_list[c(1,5,9,13)+2]
GWAS_SAS=gwas_list[c(1,5,9,13)+3]
data("AFRLDSC")
data("EASLDSC")
data("EURLDSC")
data("SASLDSC")
data("hapmap3")
GWAS_AFR=filterss(gwas_data_list=GWAS_AFR,ref_panel=hapmap3)
GWAS_SAS=filterss(gwas_data_list=GWAS_SAS,ref_panel=hapmap3)
GWAS_EUR=filterss(gwas_data_list=GWAS_EUR,ref_panel=hapmap3)
GWAS_EAS=filterss(gwas_data_list=GWAS_EAS,ref_panel=hapmap3)
fitAFR1=ldsc_poet(GWAS_List=GWAS_AFR,LDSC=AFRLDSC)
fitAFR2=ldsc_sample(GWAS_List=GWAS_AFR,LDSC=AFRLDSC)
fitEAS1=ldsc_poet(GWAS_List=GWAS_EAS,LDSC=EASLDSC)
fitEAS2=ldsc_sample(GWAS_List=GWAS_EAS,LDSC=EASLDSC)
fitEUR1=ldsc_poet(GWAS_List=GWAS_EUR,LDSC=EURLDSC)
fitEUR2=ldsc_sample(GWAS_List=GWAS_EUR,LDSC=EURLDSC)
fitSAS1=ldsc_poet(GWAS_List=GWAS_SAS,LDSC=SASLDSC)
fitSAS2=ldsc_sample(GWAS_List=GWAS_SAS,LDSC=SASLDSC)
