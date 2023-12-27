#' Estimate Heritability and Genetic Correlation Matrix Using LDSC estimated from sample LD matrix estimate
#'
#' The \code{ldscR} function estimates heritability and genetic correlation matrices using Linkage Disequilibrium Score Regression (LDSC). It processes GWAS summary statistics, harmonizes alleles with a reference panel, and computes genetic covariance and error covariance matrices.
#'
#' @param GWAS_List A list of data.frames where each data.frame contains GWAS summary statistics for a trait. Each data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
<<<<<<< HEAD
#' @param intercept.lower The lower boundary of the intercept estimate in heritability estimation.
#' @param intercept.upper The upper boundary of the intercept estimate in heritability estimation.
#' @param h2.upper The upper boundary of the heritability estimate in heritability estimation.
=======
>>>>>>> 22010f36f6745d2aa268e6c6c8c61b28157382f6
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{GCovEst}{Estimated genetic covariance matrix.}
#'   \item{GCovSE}{Standard errors of the estimated genetic covariance matrix.}
#'   \item{ECovEst}{Estimated error covariance matrix.}
#'   \item{ECovSE}{Standard errors of the estimated error covariance matrix.}
#' }
#'
#' @examples
#' data(hapmap3)
#' data(EURLDSC)
#' ref_panel <- hapmap3
#' LDSC <- EURLDSC
#' results <- ldscR(GWAS_List, ref_panel, LDSC)
#'
#' @details The \code{ldscR} function is designed for advanced genetic statistics and requires a good understanding of GWAS summary statistics, LDSC methodology, and statistical genetics. Users should ensure that input data is correctly formatted and that they understand the implications of the estimates produced by the function.
#'
#' @importFrom stats lm
#' @importFrom data.table setDT setkey copy setnames
<<<<<<< HEAD
#' @importFrom nloptr nloptr
#'
#' @export
#'
ldscR=function(GWAS_List,LDSC,intercept.lower=0.9,intercept.upper=1.1,h2.upper=0.8){
=======
#'
#' @export
#'
ldscR=function(GWAS_List,LDSC){
>>>>>>> 22010f36f6745d2aa268e6c6c8c61b28157382f6

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
  print("Processing data")
  print(t0)

  ############################# initial estimator ####################################
  t1=Sys.time()
  col_names = names(ZMatrix1)
  GCovEst1=GCovSE1=ECovEst1=ECovSE1=diag(p)*0
  for(i in 1:p){
    for(j in 1:p){
      z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
      l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[j+1]]]/M)
      fit0=lm(z~l)
      summary0=summary(fit0)
      GCovEst1[i,j]=GCovEst1[j,i]=summary0$coefficients[2,1]
      ECovEst1[i,j]=ECovEst1[j,i]=summary0$coefficients[1,1]
      GCovSE1[i,j]=GCovSE1[j,i]=summary0$coefficients[2,2]^2
      ECovSE1[i,j]=ECovSE1[j,i]=summary0$coefficients[1,2]^2
    }
  }
  t1=Sys.time()-t1
  print("Initial Genetic Covariance Estimate")
  print(t1)

  ############################ reweight for efficiency ##################################
  t2=Sys.time()
  GCovEst=GCovSE=ECovEst=ECovSE=diag(p)*0
<<<<<<< HEAD
  objective <- function(beta) {
    sum((z - X %*% beta)^2*w)
  }
  lower_bounds <- c(intercept.lower, 0)
  upper_bounds <- c(intercept.upper, h2.upper)
  for(i in 1:p){
    z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
    l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
    w=1/(1+l*GCovEst[i,i])
    fit0=lm(z~l,weights=w)
    X=cbind(1,l)
    result <- nloptr(x0 = c(1,0.1),
                     eval_f = objective,
                     lb = lower_bounds,
                     ub = upper_bounds,
                     opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
    beta=result$solution
    vare=var(c(z-X%*%beta)*sqrt(w))
    H=solve(t(X)%*%(X*sqrt(w)))*vare
    GCovEst[i,i]=beta[2]
    ECovEst[i,i]=beta[1]
    GCovSE[i,i]=sqrt(H[2,2])
    ECovSE[i,i]=sqrt(H[1,1])
=======
  for(i in 1:p){
    z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
    l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
    w=(1+l*GCovEst[i,i])
    w[w>2*median(w)]=2*median(w)
    w=w^2
    fit0=lm(z~l,weights=1/w)
    summary0=summary(fit0)
    GCovEst[i,i]=summary0$coefficients[2,1]
    ECovEst[i,i]=summary0$coefficients[1,1]
    GCovSE[i,i]=summary0$coefficients[2,2]
    ECovSE[i,i]=summary0$coefficients[1,2]
>>>>>>> 22010f36f6745d2aa268e6c6c8c61b28157382f6
  }

  for(i in 1:(p-1)){
    for(j in (i+1):p){
      z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
      l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[j+1]]]/M)
      li = LDSC1$LDSC * NMatrix1[[col_names[i+1]]]/M
      lj = LDSC1$LDSC * NMatrix1[[col_names[j+1]]]/M
      w=sqrt((1+li*GCovEst[i,i])*(1+lj*GCovEst[j,j])+(l*GCovEst1[i,j]+ECovEst1[i,j])^2)
      w[w>2*median(w)]=2*median(w)
      w=w^2
      fit0=lm(z~l,weights=1/w)
      summary0=summary(fit0)
      GCovEst[i,j]=GCovEst[j,i]=summary0$coefficients[2,1]
      ECovEst[i,j]=ECovEst[j,i]=summary0$coefficients[1,1]
      GCovSE[i,j]=GCovSE[j,i]=summary0$coefficients[2,2]
      ECovSE[i,j]=ECovSE[j,i]=summary0$coefficients[1,2]
    }
  }
  t2=Sys.time()-t2
  print("Final Genetic Covariance Estimate")
  print(t2)

  row.names(GCovEst)=colnames(GCovSE)=row.names(ECovEst)=colnames(ECovSE)=NAM
  return(A=list(GCovEst=GCovEst,GCovSE=GCovSE,ECovEst=ECovEst,ECovSE=ECovEst,stage1.time=t1,stage2.time=t2))
}
