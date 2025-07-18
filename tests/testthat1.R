data("hapmap3")
data("EURLDSC")
data("PCSK9")
PCSK9=merge_union(gwas_data_list=PCSK9,ref_panel=hapmap3,missing.thres=0.9)
Boundary=list(intercept.lower=0.5,intercept.upper=1.5,h2.upper=0.95,gcov.lower=-0.95,gcov.upper=0.95,ecov.lower=-0.95,ecov.upper=0.95)
fitPCSK9=local_ldscR(GWAS_List=PCSK9$gwas_data_list,LDSC=EURLDSC,Boundary=Boundary)
corrplot(fitPCSK9$GCovEst,method="color",is.corr=F,addCoef.col = "black",tl.cex=0.5)
corrplot(fitPCSK9$ECovEst,method="color",is.corr=F,addCoef.col = "black",tl.cex=0.5)

TIME=matrix(0,12,2)%>%as.data.frame(.)
A1=A2=matrix(0,12,4)
for(i in 1:12){
  t_1=Sys.time()
  Ajoint=merge(PCSK9$gwas_data_list[[i]],EURLDSC,by="SNP")
  ind=which(Ajoint$Zscore!=0)
  bigsnpr1=snp_ldsc(
    ld_score=Ajoint$LDSC[ind],
    ld_size=length(ind),
    chi2=Ajoint$Zscore[ind]^2,
    sample_size=Ajoint$N[ind],
    blocks = 30,
    intercept = NULL,
    chi2_thr1 = 30,
    chi2_thr2 = Inf,
    ncores = 1
  )
  t_1=Sys.time()-t_1

  t_2=Sys.time()
  bigsnpr2=snp_ldsc(
    ld_score=Ajoint$LDSC[ind],
    ld_size=length(Ajoint$SNP[ind]),
    chi2=Ajoint$Zscore[ind]^2,
    sample_size=Ajoint$N[ind],
    blocks = 50,
    intercept = 1,
    chi2_thr1 = 30,
    chi2_thr2 = Inf,
    ncores = 1
  )
  t_2=Sys.time()-t_2
  A1[i,]=bigsnpr1
  A2[i,]=bigsnpr2
  TIME[i,]=c(t_1,t_2)
  print(TIME[i,])
}
