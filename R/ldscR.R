#' Estimate Heritability and Genetic Correlation Matrix Using LDSC estimated from sample LD matrix estimate
#'
#' The \code{ldscR} function estimates heritability and genetic correlation matrices using Linkage Disequilibrium Score Regression (LDSC). It processes GWAS summary statistics, harmonizes alleles with a reference panel, and computes genetic covariance and error covariance matrices.
#'
#' @param GWAS_List A list of data.frames where each data.frame contains GWAS summary statistics for a trait. Each data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param Boundary A list of the information of upper and lower boundaries of parameter estimates. It typically includes intercept.lower (he lower boundary of the intercept estimate), intercept.upper (the upper boundary of the intercept estimate), h2.upper (the upper boundary of the heritability estimate).
#' @param zsquare_thresh Threshold on the z-score square in heritability estimation.
#' @param cov_thresh Threshold on the abs(z-score1 * z-score2) in genetic covariance estimation.
#' @param estimate_SE If estimating the standard error of the genetic covariance matrix and estimation error covariance matrix. Default to F.
#' @param nblock The number of blocks for bootstrap-based standard error estimation. Default to 200.
#' @param sampling.time time for the estimation of standard errors. Default to 200.
#' @param sampling.ratio The sub-sampling ratio in each bootstrap. Default to 0.5.

#' @return A list containing the following elements:
#' \itemize{
#'   \item{GCovEst}{Estimated genetic covariance matrix.}
#'   \item{GCovSE}{Standard errors of the estimated genetic covariance matrix.}
#'   \item{ECovEst}{Estimated error covariance matrix.}
#'   \item{ECovSE}{Standard errors of the estimated error covariance matrix.}
#'   \item{Computing.time}{Computing time in each stage.}
#' }
#'
#' @examples
#' data(hapmap3)
#' data(EURLDSC)
#' ref_panel <- hapmap3
#' LDSC <- EURLDSC
#' results <- ldscR(GWAS_List,LDSC)
#'
#' @details The \code{ldscR} function is designed for advanced genetic statistics and requires a good understanding of GWAS summary statistics, LDSC methodology, and statistical genetics. Users should ensure that input data is correctly formatted and that they understand the implications of the estimates produced by the function.
#'
#' @importFrom stats lm cov2cor runif sd
#' @import data.table
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom nloptr nloptr
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#'
#'
ldscR=function(GWAS_List,LDSC,Boundary=F,zsquare_thresh=50000,cov_thresh=10000,estimate_SE=F,nblock=500,sampling.time=500,sampling.ratio=0.5){

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
if(Boundary[1]==F){
Boundary=list(intercept.lower=0.95,intercept.upper=1.05,h2.upper=0.95)
}
t0 = Sys.time() - t0
cat("Processing data -> ")
print(t0)

############################# initial estimator ####################################
t1=Sys.time()
col_names = names(ZMatrix1)
GCovEst1 = ECovEst1 = diag(p)*0
objective <- function(beta) {
sum((z - X %*% beta)^2)
}
lower_bounds <- c(Boundary$intercept.lower, 0)
upper_bounds <- c(Boundary$intercept.upper, Boundary$h2.upper)

for(i in 1:p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
z[which(z>zsquare_thresh)] = zsquare_thresh
X=cbind(1,l)
result <- nloptr(x0 = c(1,0.1),
eval_f = objective,
lb = lower_bounds,
ub = upper_bounds,
opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
beta=result$solution
GCovEst1[i,i]=beta[2]
ECovEst1[i,i]=beta[1]
}

for(i in 1:(p-1)){
for(j in (i+1):p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
z[which(abs(z)>cov_thresh)] = sign(z[which(abs(z)>cov_thresh)]) * cov_thresh
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[j+1]]]/M)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X)
Xty=matrixVectorMultiply(t(X),z)
beta=c(solve(XtX)%*%Xty)
GCovEst1[i,j]=GCovEst1[j,i]=beta[2]
ECovEst1[i,j]=ECovEst1[j,i]=beta[1]
}
}
t1=Sys.time()-t1
cat("Initial Genetic Covariance Estimate -> ")
print(t1)

############################ reweight for efficiency ##################################
t2=Sys.time()
objective <- function(beta) {
sum((z - X %*% beta)^2*w)
}
GCovEst=ECovEst=diag(p)*0
lower_bounds <- c(Boundary$intercept.lower, 0)
upper_bounds <- c(Boundary$intercept.upper, Boundary$h2.upper)

for(i in 1:p){
z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
z[which(z>zsquare_thresh)] = zsquare_thresh
l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
w = 1/(1+l*GCovEst[i,i])^2
X=cbind(1,l)
result <- nloptr(x0 = c(1,0.1),
eval_f = objective,
lb = lower_bounds,
ub = upper_bounds,
opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
beta=result$solution
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
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=c(solve(XtX)%*%Xty)
GCovEst[i,j]=GCovEst[j,i]=beta[2]
ECovEst[i,j]=ECovEst[j,i]=beta[1]
}
}
t2=Sys.time()-t2
cat("Final Genetic Covariance Estimate -> ")
print(t2)


############################ resampling for standard error ##################################
if(estimate_SE==T){
t3=Sys.time()
objective <- function(beta) {
sum((z - X %*% beta)^2*w)
}
lower_bounds <- c(Boundary$intercept.lower, 0)
upper_bounds <- c(Boundary$intercept.upper, Boundary$h2.upper)

blocksize=floor(M/nblock)
IndexMatrix=matrix(c(1:(nblock*blocksize)),blocksize,nblock)
gap=M-(nblock*blocksize)
GCovEstt=ECovEstt=GCorEstt=array(0,c(sampling.time,p,p))
GCovSE=ECovSE=GCorSE=diag(p)*0
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(t in 1:sampling.time){
setTxtProgressBar(pb, t)
indcol=sample(ncol(IndexMatrix),ncol(IndexMatrix)*sampling.ratio,replace=F)
ind=as.vector(IndexMatrix[,indcol])+floor(runif(1,1,gap))
ZMatrixt=ZMatrix1[ind,]
NMatrixt=NMatrix1[ind,]
LDSCt=LDSC1[ind,]
for(i in 1:p){
z = ZMatrixt[[col_names[i+1]]] * ZMatrixt[[col_names[i+1]]]
z[which(z>zsquare_thresh)] = zsquare_thresh
l = LDSCt$LDSC * sqrt(NMatrixt[[col_names[i+1]]]/M) * sqrt(NMatrixt[[col_names[i+1]]]/M)
w = 1/(1+l*GCovEst[i,i])^2
X=cbind(1,l)
result <- nloptr(x0 = c(1,0.1),
eval_f = objective,
lb = lower_bounds,
ub = upper_bounds,
opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
beta=result$solution
GCovEstt[t,i,i]=beta[2]
ECovEstt[t,i,i]=beta[1]
}

for(i in 1:(p-1)){
for(j in (i+1):p){
z = ZMatrixt[[col_names[i+1]]] * ZMatrixt[[col_names[j+1]]]
z[which(abs(z)>cov_thresh)] = sign(z[which(abs(z)>cov_thresh)]) * cov_thresh
l = LDSCt$LDSC * sqrt(NMatrixt[[col_names[i+1]]]/M) * sqrt(NMatrixt[[col_names[j+1]]]/M)
li = LDSCt$LDSC * NMatrixt[[col_names[i+1]]]/M
lj = LDSCt$LDSC * NMatrixt[[col_names[j+1]]]/M
w = 1/((1+li*GCovEst1[i,i])*(1+lj*GCovEst1[j,j])+(l*GCovEst1[i,j]+ECovEst1[i,j])^2)
X=cbind(1,l)
XtX=matrixMultiply(t(X),X*w)
Xty=matrixVectorMultiply(t(X),z*w)
beta=c(solve(XtX)%*%Xty)
GCovEstt[t,i,j]=GCovEstt[t,j,i]=beta[2]
ECovEstt[t,i,j]=ECovEstt[t,j,i]=beta[1]
}
}
GCorEstt[t,,]=cov2cor(GCovEstt[t,,])
}
close(pb)
for(i in 1:p){
for(j in i:p){
GCovSE[i,j]=GCovSE[j,i]=sd(GCovEstt[,i,j])
ECovSE[i,j]=ECovSE[j,i]=sd(ECovEstt[,i,j])
GCorSE[i,j]=GCorSE[j,i]=sd(GCorEstt[,i,j])
}
}
t3=Sys.time()-t3
cat("Resampling for Standard Error -> ")
print(t3)
}else{
GCovSE=ECovSE=GCorSE=0*diag(p)
t3=0
}
GCorEst=cov2cor(GCovEst)
GCorEst[is.na(GCorEst)]=0
row.names(GCovEst)=colnames(GCovEst)=row.names(GCorEst)=colnames(GCorEst)=row.names(ECovEst)=colnames(ECovEst)=NAM
row.names(GCovSE)=colnames(GCovSE)=row.names(GCorSE)=colnames(GCorSE)=row.names(ECovSE)=colnames(ECovSE)=NAM
row.names(GCovEst1)=colnames(GCovEst1)=row.names(ECovEst1)=colnames(ECovEst1)=NAM
Computing.time=list(data_organization=t0,initial_estimation=t1,reweighting_estimation=t2,resampling_estimation=t3)
return(A=list(GCovEst=GCovEst,GCovSE=GCovSE,GCorEst=GCorEst,GCorSE=GCorSE,ECovEst=ECovEst,ECovSE=ECovSE,Computing.time=Computing.time))
}

