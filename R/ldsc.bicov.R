#' Single-Variate Linkage Disequilibrium Score Regression (LDSC)
#'
#' The `ldsc.bicov` function performs single-variate Linkage Disequilibrium Score Regression (LDSC) analysis. It is designed to estimate genetic covariance from two GWAS summary statistics, accounting for linkage disequilibrium (LD) between SNPs.
#'
#' @param gwas1 A data.frame containing GWAS summary statistics for the first trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param gwas2 A data.frame containing GWAS summary statistics for the second trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param h1 The heritability of the first trait.
#' @param h2 The heritability of the second trait.
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param zsquare_thresh A threshold for the squared Z-scores in heritability estimation to control for extreme values.
#' @param nblock The number of blocks for bootstrap-based standard error estimation.
#' @param sampling.time The number of block bootstrap.
#' @param sampling.ratio The sub-sampling ratio in each bootstrap.
#' @return A data.frame containing heritability estimates and their standard errors, along with the intercept and its standard error.
#'
#' @importFrom stats lm
#' @importFrom data.table setDT setkey copy setnames
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#'
#' @export
#'
ldsc.bicov=function(gwas1,gwas2,h1,h2,LDSC,zsquare_thresh=50000,nblock=500,sampling.time=200,sampling.ratio=0.5){

############################# Basic Information ###############################

t0 = Sys.time()
LDSC$Order=1:nrow(LDSC)
gwas1=merge(gwas1,LDSC,by="SNP")
gwas2=merge(gwas2,LDSC,by="SNP")
snp_intersect=intersect(gwas1$SNP,gwas2$SNP)
gwas1=gwas1[which(gwas1$SNP%in%snp_intersect),]
gwas2=gwas2[which(gwas2$SNP%in%snp_intersect),]
gwas1 <- gwas1[order(gwas1$Order), ]
gwas2 <- gwas2[order(gwas2$Order), ]
M=nrow(gwas1)
t0 = Sys.time() - t0
cat("Processing data -> ")
print(t0)

############################# initial estimator ####################################
t1=Sys.time()
z = gwas1$Zscore*gwas2$Zscore
l = gwas1$LDSC/M*sqrt(gwas1$N)*sqrt(gwas2$N)
z[which(abs(z)>zsquare_thresh)] = zsquare_thresh*sign(z[which(abs(z)>zsquare_thresh)])
X=cbind(1,l)
XtX=matrixMultiply(t(X),X)
Xty=matrixVectorMultiply(t(X),z)
theta=c(solve(XtX)%*%Xty)
t1=Sys.time()-t1
cat("Initial Covariance Estimate -> ")
print(t1)

############################ reweighting for efficiency ##################################
t2=Sys.time()
z = gwas1$Zscore*gwas2$Zscore
l = gwas1$LDSC/M*sqrt(gwas1$N)*sqrt(gwas2$N)
z[which(abs(z)>zsquare_thresh)] = zsquare_thresh*sign(z[which(abs(z)>zsquare_thresh)])
X=cbind(1,l)
w=c((gwas1$LDSC/M*gwas1$N)*(gwas2$LDSC/M*gwas2$N)+(l*theta[2]+theta[1])^2)
w=1/w
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
theta=c(solve(XtX)%*%Xty)
t2=Sys.time()-t2
cat("Genetic Covariance Estimation -> ")
print(t2)

########################### resampling for standard error ############################
t3=Sys.time()
w=c((gwas1$LDSC/M*gwas1$N)*(gwas2$LDSC/M*gwas2$N)+(l*theta[2]+theta[1])^2)
w=1/w
blocksize=floor(M/nblock)
IndexMatrix=matrix(c(1:(nblock*blocksize)),blocksize,nblock)
gap=M-(nblock*blocksize)
h2.vec=intercept.vec=c(1:sampling.time)
for(i in 1:sampling.time){
indcol=sample(ncol(IndexMatrix),ncol(IndexMatrix)*sampling.ratio,replace=F)
ind=as.vector(IndexMatrix[,indcol])+floor(runif(1,1,gap))
z1=z[ind]
l1=l[ind]
w1=w[ind]
X1=cbind(1,l1)
XtX1=matrixMultiply(t(X1),X1*w1)
Xty1=matrixVectorMultiply(t(X1),z1*w1)
theta1=c(solve(XtX1)%*%Xty1)
h2.vec[i]=theta1[2]
intercept.vec[i]=theta1[1]
}
t3=Sys.time()-t3
cat("Standard Error Estimation -> ")
print(t3)
ecov.se=sqrt(mean((intercept.vec-theta[1])^2))
gcov.se=sqrt(mean((h2.vec-theta[2])^2))
return(A=data.frame(ecov=theta[1],ecov.se=ecov.se,gcov=theta[2],gcov.se=gcov.se))
}

