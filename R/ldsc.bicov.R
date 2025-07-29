#' Single-Variate Linkage Disequilibrium Score Regression (LDSC)
#'
#' The `ldsc.bicov` function performs bi-variate genetic covariance using Linkage Disequilibrium Score Regression (LDSC) analysis. It is designed to estimate heritability from gwas1 summary statistics, accounting for linkage disequilibrium (LD) between SNPs. The function harmonizes gwas1 data with LD scores and applies non-linear optimization to estimate heritability.
#'
#' @param gwas1 A data.frame containing gwas summary statistics for the first trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N).
#' @param gwas2 A data.frame containing gwas summary statistics for the second trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param h21 The heritability estimate of the first trait.
#' @param h22 The heritability estimate of the second trait.
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param nblock The number of blocks for bootstrap-based standard error estimation.
#' @param sampling.time The number of block bootstrap. If sampling.time=0, only genetic covariance estimate will be returned.
#' @return A data.frame containing genetic covariance estimates and their standard errors, along with the estimation error covariance estimate and its standard error.

#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom data.table setDT setkey copy setnames
#'
#' @export
#'
ldsc.bicov=function(gwas1,gwas2,h21,h22,LDSC,nblock=500,sampling.time=0){
setDT(LDSC)
setkey(LDSC, SNP)

############################# Basic Information ###############################
t0 = Sys.time()
SNPintersect=intersect(gwas1$SNP,gwas2$SNP)
gwas1=gwas1[which(gwas1$SNP%in%SNPintersect),]
gwas2=gwas2[which(gwas2$SNP%in%SNPintersect),]
gwas1 = LDSC[gwas1, nomatch=0]
gwas2 = LDSC[gwas2, nomatch=0]
M=nrow(gwas1)
t0 = Sys.time() - t0
print("Processing data")
print(t0)

############################# initial estimator ####################################
t1=Sys.time()
z = gwas1$Zscore*gwas2$Zscore
l = gwas1$LDSC*sqrt(gwas1$N)*sqrt(gwas2$N)/M
w = 1/(1+gwas1$LDSC*gwas1$N/M*h21)/(1+gwas2$LDSC*gwas2$N/M*h22)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=as.vector(solve(XtX)%*%Xty)
gcov.ini=beta[2]
ecov.ini=beta[1]
t1=Sys.time()-t1
print("Initial Genetic Covariance Estimate")
print(t1)

############################ reweight for efficiency ##################################
t2=Sys.time()
z = gwas1$Zscore*gwas2$Zscore
l = gwas1$LDSC*sqrt(gwas1$N)*sqrt(gwas2$N)/M
w = (1+gwas1$LDSC*gwas1$N/M*h21)/(1+gwas2$LDSC*gwas2$N/M*h22)+2*(l*gcov.ini+ecov.ini)^2
w = 1/w
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=as.vector(solve(XtX)%*%Xty)
gcov=beta[2]
ecov=beta[1]

w = (1+gwas1$LDSC*gwas1$N/M*h21)/(1+gwas2$LDSC*gwas2$N/M*h22)+2*(l*gcov+ecov)^2
w = 1/w
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=as.vector(solve(XtX)%*%Xty)
gcov=beta[2]
ecov=beta[1]
t2=Sys.time()-t2
print("Heritability Estimation")
print(t2)

################################################################################
if(sampling.time>0){
t3=Sys.time()
blocksize=floor(M/nblock)
IndexMatrix=matrix(c(1:(nblock*blocksize)),blocksize,nblock)
gap=M-(nblock*blocksize)
h2.vec=intercept.vec=c(1:sampling.time)
for(i in 1:sampling.time){
indcol=sample(ncol(IndexMatrix),ncol(IndexMatrix),replace=T)
ind=as.vector(IndexMatrix[,indcol])+ifelse(gap>1,floor(runif(1,0,gap)),0)
z1=z[ind]
l1=l[ind]
w1=w[ind]
X1=cbind(1,l1)
XtX=matrixMultiply(t(X1),X1*w1)
Xty=matrixVectorMultiply(t(X1),z1*w1)
beta=as.vector(solve(XtX)%*%Xty)
h2.vec[i]=beta[2]
intercept.vec[i]=beta[1]
}
t3=Sys.time()-t3
print("Standard Error Estimation")
print(t3)
ecov.se=sqrt(mean((intercept.vec-ecov)^2))
gcov.se=sqrt(mean((h2.vec-gcov)^2))
}else{
gcov.se=ecov.se=NA
}
return(A=data.frame(ecov=ecov,ecov.se=ecov.se,gcov=gcov,gcov.se=gcov.se,M=M))
}

