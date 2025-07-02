#' Estimate Heritability and Genetic Correlation Matrix Using LDSC estimated from sample LD matrix estimate
#'
#' The \code{ldscR} function estimates heritability and genetic correlation matrices using Linkage Disequilibrium Score Regression (LDSC). It processes GWAS summary statistics, harmonizes alleles with a reference panel, and computes genetic covariance and error covariance matrices.
#'
#' @param GWAS_List A list of data.frames where each data.frame contains GWAS summary statistics for a trait. Each data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param zsquare_thresh Threshold on the z-score square in heritability estimation.
#' @param cov_thresh Threshold on the abs(z-score1 * z-score2) in genetic covariance estimation.
#' @param estimate_SE If estimating the standard error of the genetic covariance matrix and estimation error covariance matrix. Default to F.
#' @param sampling.time Bootstrapping time for the estimation of standard errors. Default to 200.
#' @param n_threads Number of threads used in Block-wise bootstrap. Default to NULL which uses half of CPU cores detected via OpenMP.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{GCovEst}{Estimated genetic covariance matrix.}
#'   \item{GCovSE}{Standard errors of the estimated genetic covariance matrix.}
#'   \item{ECovEst}{Estimated error covariance matrix.}
#'   \item{ECovSE}{Standard errors of the estimated error covariance matrix.}
#'   \item{Computing.time}{Computing time in each stage.}
#' }
#'
#' @importFrom stats lm cov2cor runif sd
#' @import data.table
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#'
#' @export
#'
#'
ldscR=function(GWAS_List,LDSC,zsquare_thresh=50000,cov_thresh=10000,estimate_SE=F,sampling.time=500,n_threads=NULL){

############################# Basic Information ###############################
t0 = Sys.time()
GWAS_List_dt = lapply(GWAS_List, setDT)
NAM = names(GWAS_List)
p = length(GWAS_List)
ZMatrix = copy(GWAS_List_dt[[1]])[, .(SNP, Zscore = Zscore)]
setnames(ZMatrix, "Zscore", NAM[1])
NMatrix = copy(GWAS_List_dt[[1]])[, .(SNP, N = N)]
setnames(NMatrix, "N", NAM[1])
for (i in 2:p) {
ZMatrix[, (NAM[i]) := GWAS_List_dt[[i]][ZMatrix, on = .(SNP), x.Zscore]]
NMatrix[, (NAM[i]) := GWAS_List_dt[[i]][NMatrix, on = .(SNP), x.N]]
}
rm(GWAS_List_dt)
setkey(ZMatrix, SNP)
setkey(NMatrix, SNP)
LDSC_dt = setDT(LDSC)
snplist = intersect(ZMatrix$SNP, LDSC_dt$SNP)
ZMatrix1 = ZMatrix[snplist,]
NMatrix1 = NMatrix[snplist,]
LDSC1 = LDSC_dt[snplist, on = .(SNP)]
M = length(snplist)
rm(ZMatrix, NMatrix, LDSC_dt, snplist)
t0 = Sys.time() - t0
cat("Processing data -> ")
print(t0)

GCovEst1=GCovEst=ECovEst=ECovEst1=diag(p)
############################# initial estimator ####################################
t1=Sys.time()
col_names = names(ZMatrix1)
for(i in 1:p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
z[which(z>zsquare_thresh)] = zsquare_thresh
X=cbind(1,l)
XtX=matrixMultiply(t(X),X)/nrow(X)
Xty=matrixVectorMultiply(t(X),z)/nrow(X)
beta=matrixVectorMultiply(solve(XtX),Xty)
GCovEst1[i,i]=beta[2]
ECovEst1[i,i]=beta[1]
}

for(i in 1:(p-1)){
for(j in (i+1):p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
z[which(abs(z)>cov_thresh)] = sign(z[which(abs(z)>cov_thresh)]) * cov_thresh
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[j+1]]]/M)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X)/nrow(X)
Xty=matrixVectorMultiply(t(X),z)/nrow(X)
beta=matrixVectorMultiply(solve(XtX),Xty)
GCovEst1[i,j]=GCovEst1[j,i]=beta[2]
ECovEst1[i,j]=ECovEst1[j,i]=beta[1]
}
}
t1=Sys.time()-t1
cat("Initial Genetic Covariance Estimate -> ")
print(t1)

############################ reweight for efficiency ##################################
t2=Sys.time()

for(i in 1:p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
z[which(z>zsquare_thresh)] = zsquare_thresh
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
w = 1/(1+l*GCovEst[i,i])^2
XtX=matrixMultiply(t(X),X*w)/nrow(X)
Xty=matrixVectorMultiply(t(X),z*w)/nrow(X)
beta=matrixVectorMultiply(solve(XtX),Xty)
GCovEst[i,i]=beta[2]
ECovEst[i,i]=beta[1]
}

for(i in 1:(p-1)){
for(j in (i+1):p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
z[which(abs(z)>cov_thresh)] = sign(z[which(abs(z)>cov_thresh)]) * cov_thresh
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[j+1]]]/M)
li = LDSC1$LDSC * NMatrix1[[col_names[i+1]]]/M
lj = LDSC1$LDSC * NMatrix1[[col_names[j+1]]]/M
w = 1/((1+li*GCovEst1[i,i])*(1+lj*GCovEst1[j,j])+(l*GCovEst1[i,j]+ECovEst1[i,j])^2)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)/nrow(X)
Xty=matrixVectorMultiply(t(X),z*w)/nrow(X)
beta=matrixVectorMultiply(solve(XtX),Xty)
GCovEst[i,j]=GCovEst[j,i]=beta[2]
ECovEst[i,j]=ECovEst[j,i]=beta[1]
}
}
t2=Sys.time()-t2
cat("Final Genetic Covariance Estimate -> ")
print(t2)

############################ resampling for standard error ##################################
t3=Sys.time()
if(estimate_SE==T){
bootstrap_res <- ldsc_bootstrap_wrapper(ZMatrix1 = ZMatrix1,
                                      NMatrix1 = NMatrix1,
                                      LDSC1 = LDSC1,
                                      GCovInit = GCovEst,
                                      ECovInit = ECovEst,
                                      sampling.time = sampling.time,
                                      n_threads=n_threads)
GCovSE <- bootstrap_res$GCovSE
ECovSE <- bootstrap_res$ECovSE
GCorSE <- bootstrap_res$GCorSE
}else{
GCovSE=ECovSE=GCorSE=diag(p)*0
}
t3=Sys.time()-t3

GCorEst=cov2cor_ldscr(GCovEst)
row.names(GCovEst)=colnames(GCovEst)=row.names(GCorEst)=colnames(GCorEst)=row.names(ECovEst)=colnames(ECovEst)=NAM
row.names(GCovSE)=colnames(GCovSE)=row.names(GCorSE)=colnames(GCorSE)=row.names(ECovSE)=colnames(ECovSE)=NAM
row.names(GCovEst1)=colnames(GCovEst1)=row.names(ECovEst1)=colnames(ECovEst1)=NAM
Computing.time=list(data_organization=t0,initial_estimation=t1,reweighting_estimation=t2,resampling_estimation=t3)
return(A=list(GCovEst=GCovEst,GCovSE=GCovSE,GCorEst=GCorEst,GCorSE=GCorSE,ECovEst=ECovEst,ECovSE=ECovSE,Computing.time=Computing.time))
}

