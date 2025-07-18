#' Estimate Heritability and Genetic Correlation Matrix Using LDSC estimated from sample LD matrix estimate
#'
#' The \code{ldscR} function estimates heritability and genetic correlation matrices using Linkage Disequilibrium Score Regression (LDSC). It processes GWAS summary statistics, harmonizes alleles with a reference panel, and computes genetic covariance and error covariance matrices.
#'
#' @param GWAS_List A list of data.frames where each data.frame contains GWAS summary statistics for a trait. Each data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param Boundary A list of the information of upper and lower boundaries of parameter estimates. It typically includes intercept.lower (he lower boundary of the intercept estimate), intercept.upper (the upper boundary of the intercept estimate), h2.upper (the upper boundary of the heritability estimate), ecov.lower (the lower boundary of the estimation error covariance estimate), ecov.upper (the upper boundary of the estimation error covariance estimate), gcov.lower (the lower boundary of the genetic covariance estimate), gcov.upper (the upper boundary of the genetic covariance estimate).
#' @param zsquare_thresh Threshold on the z-score square in heritability estimation.
#' @param cov_thresh Threshold on the abs(z-score1 * z-score2) in genetic covariance estimation.
#' @param min.eps  A lower boundary of the eigenvalues in genetic correlation matrix estimate and estimation error correlation matrix estimate.
#' @return A list containing the following elements:
#' \itemize{
#'   \item{GCovEst}{Estimated genetic covariance matrix.}
#'   \item{GCovSE}{Standard errors of the estimated genetic covariance matrix.}
#'   \item{ECovEst}{Estimated error covariance matrix.}
#'   \item{ECovSE}{Standard errors of the estimated error covariance matrix.}
#'   \item{Estimate.ini}{Initial estimates with all weights being 1.}
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
#' @importFrom stats lm
#' @importFrom data.table setDT setkey copy setnames
#' @importFrom nloptr nloptr
#'
#' @export
#'
#'
ldscR=function(GWAS_List,LDSC,Boundary=F,zsquare_thresh=500,cov_thresh=300){

############################# Basic Information ###############################
t0 = Sys.time()
GWAS_List_dt = lapply(GWAS_List, setDT)
NAM = names(GWAS_List)
p = length(GWAS_List)
SNPOrder = data.table(SNP = GWAS_List_dt[[1]]$SNP, Order = 1:nrow(GWAS_List_dt[[1]]))
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
ZMatrix1 = ZMatrix1[SNPOrder[ZMatrix1, on = .(SNP), nomatch = 0]$Order]
NMatrix1 = NMatrix1[SNPOrder[NMatrix1, on = .(SNP), nomatch = 0]$Order]
LDSC1 = LDSC1[SNPOrder[LDSC1, on = .(SNP), nomatch = 0]$Order]
M = length(snplist)
rm(ZMatrix, NMatrix, LDSC_dt, snplist, SNPOrder)
t0 = Sys.time() - t0
cat("Processing data")
print(t0)

############################# initial estimator ####################################
t1=Sys.time()

ZMatrix2=as.matrix(ZMatrix1[,-1])
NMatrix2=as.matrix(NMatrix1[,-1])
NMatrix2=vectorMatrixDotProduct(LDSC1$LDSC/M,NMatrix2)
Zself=selfPairwiseProduct(ZMatrix2)
Nself=selfPairwiseProduct(NMatrix2)

fit1=ldsc(selfProductZ=Zself$selfProducts,pairwiseProductZ=Zself$pairwiseProducts,selfProductN=Nself$selfProducts,pairwiseProductN=Nself$pairwiseProducts)
fit1$GCovEst1=fit1$GCovEst1+t(fit1$GCovEst1)-diag(diag(fit1$GCovEst1))
fit1$ECovEst1=fit1$ECovEst1+t(fit1$ECovEst1)-diag(diag(fit1$ECovEst1))
t1=Sys.time()-t1
cat("Initial Genetic Covariance Estimate")
print(t1)

############################ reweight for efficiency ##################################
t2=Sys.time()

t2=Sys.time()-t2
print("Final Genetic Covariance Estimate")
print(t2)

row.names(GCovEst)=colnames(GCovEst)=row.names(ECovEst)=colnames(ECovEst)=NAM
row.names(GCovSE)=colnames(GCovSE)=row.names(ECovSE)=colnames(ECovSE)=NAM
row.names(GCovEst1)=colnames(GCovEst1)=row.names(ECovEst1)=colnames(ECovEst1)=NAM
Estimate.ini=list(GCovEst=GCovEst1,ECovEst=ECovEst1)
Computing.time=list(stage1.time=t0,stage2.time=t1,stage3.time=t2)
return(A=list(GCovEst=GCovEst,GCovSE=GCovSE,ECovEst=ECovEst,ECovSE=ECovSE,Estimate.ini=Estimate.ini,Computing.time=Computing.time))
}

