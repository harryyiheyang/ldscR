#' Single-Variate Linkage Disequilibrium Score Regression (LDSC)
#'
#' The `ldsc.uni` function performs single-variate Linkage Disequilibrium Score Regression (LDSC) analysis. It is designed to estimate heritability from GWAS summary statistics, accounting for linkage disequilibrium (LD) between SNPs. The function harmonizes GWAS data with LD scores and applies non-linear optimization to estimate heritability.
#'
#' @param gwas A data.frame containing GWAS summary statistics for a single trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param nblock The number of blocks for bootstrap-based standard error estimation.
#' @param sampling.time The number of block bootstrap. If sampling.time=0, only genetic covariance estimate will be returned.
#' @return A data.frame containing heritability estimates and their standard errors, along with the intercept and its standard error.

#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom data.table setDT setkey copy setnames
#'
#' @export
#'
ldsc.univ=function(gwas,LDSC,nblock=500,sampling.time=500){

############################# Basic Information ###############################
t0 = Sys.time()
LDSC$Order=1:nrow(LDSC)
gwas=merge(gwas,LDSC,by="SNP")
gwas <- gwas[order(gwas$Order), ]
M=nrow(gwas)
t0 = Sys.time() - t0
print("Processing data")
print(t0)

############################# initial estimator ####################################
t1=Sys.time()
z = gwas$Zscore^2
l = gwas$LDSC*gwas$N/M
w = 1/(1+l*0.1)^2
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=as.vector(solve(XtX)%*%Xty)
h2.ini=beta[2]
intercept.ini=beta[1]
t1=Sys.time()-t1
print("Initial Genetic Covariance Estimate")
print(t1)

############################ reweight for efficiency ##################################
t2=Sys.time()
z = gwas$Zscore^2
l = gwas$LDSC*gwas$N/M
w = 1/(1+l*h2.ini)^2
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=as.vector(solve(XtX)%*%Xty)
h2=beta[2]
intercept=beta[1]

w = 1/(1+l*h2)^2
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=as.vector(solve(XtX)%*%Xty)
h2=beta[2]
intercept=beta[1]
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
intercept.se=sqrt(mean((intercept.vec-intercept)^2))
h2.se=sqrt(mean((h2.vec-h2)^2))
}else{
h2.se=intercept.se=NA
}
return(A=data.frame(intercept=intercept,intercept.se=intercept.se,h2=h2,h2.se=h2.se,M=M))
}

