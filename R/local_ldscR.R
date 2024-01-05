#' Estimate Heritability and Genetic Correlation Matrix Using LDSC estimated from sample LD matrix estimate in a genome region.
#'
#' The \code{local_ldscR} function estimates heritability and genetic correlation matrices using Linkage Disequilibrium Score Regression (LDSC). It processes GWAS summary statistics, harmonizes alleles with a reference panel, and computes genetic covariance and error covariance matrices.
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
#' results <- local_ldscR(GWAS_List,LDSC)
#'
#' @details The \code{local_ldscR} function is designed for advanced genetic statistics and requires a good understanding of GWAS summary statistics, LDSC methodology, and statistical genetics. Users should ensure that input data is correctly formatted and that they understand the implications of the estimates produced by the function. The definition of genetic correlation is from Zhao and Zhu (2022).
#'
#' @importFrom stats lm
#' @importFrom data.table setDT setkey copy setnames
#' @importFrom nloptr nloptr
#'
#' @export
#'
#' @references Zhao B, Zhu H. On genetic correlation estimation with summary statistics from genome-wide association studies. Journal of the American Statistical Association. 2022 Jan 2;117(537):1-1.
#'
local_ldscR=function(GWAS_List,LDSC,Boundary=F,zsquare_thresh=50,cov_thresh=50,min.eps=0.001){
  
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
    Boundary=list(intercept.lower=0.95,intercept.upper=1.05,h2.lower=0.005,h2.upper=0.95,gcov.lower=-0.95,gcov.upper=0.95,ecov.lower=-0.95,ecov.upper=0.95)
  }
  t0 = Sys.time() - t0
  print("Processing data")
  print(t0)
  
  ############################# initial estimator ####################################
  t1=Sys.time()
  col_names = names(ZMatrix1)
  GCovEst1 = GCovSE1 = ECovEst1 = ECovSE1 = diag(p)*0
  objective <- function(beta) {
    sum((z - X %*% beta)^2)
  }
  lower_bounds <- c(Boundary$intercept.lower, Boundary$h2.lower)
  upper_bounds <- c(Boundary$intercept.upper, Boundary$h2.upper)
  
  for(i in 1:p){
    ind = which(ZMatrix1[[col_names[i+1]]]!=0)
    M1 = length(ind)
    z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
    l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M1) * sqrt(NMatrix1[[col_names[i+1]]]/M1)
    l = l[ind]
    z = z[ind]
    z[which(z>zsquare_thresh)] = zsquare_thresh
    X=cbind(1,l)
    result <- nloptr(x0 = c(1,0.1),
                     eval_f = objective,
                     lb = lower_bounds,
                     ub = upper_bounds,
                     opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
    beta=result$solution
    vare=var(c(z-X%*%beta))
    H=solve(t(X)%*%X)*vare
    GCovEst1[i,i]=beta[2]
    ECovEst1[i,i]=beta[1]
    GCovSE1[i,i]=sqrt(H[2,2])
    ECovSE1[i,i]=sqrt(H[1,1])
  }
  
  lower_bounds <- c(Boundary$ecov.lower, Boundary$gcov.lower)
  upper_bounds <- c(Boundary$ecov.upper, Boundary$gcov.upper)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      ind = union(which(ZMatrix1[[col_names[i+1]]]!=0),which(ZMatrix1[[col_names[j+1]]]!=0))
      M1 = length(M1)
      z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
      l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M1) * sqrt(NMatrix1[[col_names[j+1]]]/M1)
      z = z[ind]
      l = l[ind]
      z[which(abs(z)>cov_thresh)] = sign(z[which(abs(z)>cov_thresh)]) * cov_thresh
      X=cbind(1,l)
      a=c(sqrt(ECovEst1[i,i]*ECovEst1[j,j]),sqrt(GCovEst1[i,i]*GCovEst1[j,j]))
      result <- nloptr(x0 = c(0,0),
                       eval_f = objective,
                       lb = lower_bounds*a,
                       ub = upper_bounds*a,
                       opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
      beta=result$solution
      vare=var(c(z-X%*%beta))
      H=solve(t(X)%*%X)*vare
      GCovEst1[i,j]=GCovEst1[j,i]=beta[2]
      ECovEst1[i,j]=ECovEst1[j,i]=beta[1]
      GCovSE1[i,j]=GCovSE1[j,i]=sqrt(H[2,2])
      ECovSE1[i,j]=ECovSE1[j,i]=sqrt(H[1,1])
    }
  }
  t1=Sys.time()-t1
  print("Initial Genetic Covariance Estimate")
  print(t1)
  
  ############################ reweight for efficiency ##################################
  t2=Sys.time()
  
  objective <- function(beta) {
    sum((z - X %*% beta)^2*w)
  }
  GCovEst = GCovSE = ECovEst = ECovSE = diag(p)*0
  ECovEst11 <- positive_adj(cov2cor(ECovEst1),min.eps=min.eps)
  GCovEst11 <- positive_adj(GCovEst1,min.eps=min.eps)
  lower_bounds <- c(Boundary$intercept.lower, Boundary$h2.lower)
  upper_bounds <- c(Boundary$intercept.upper, Boundary$h2.upper)
  
  for(i in 1:p){
    ind = which(ZMatrix1[[col_names[i+1]]]!=0)
    M1 = length(ind)
    z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[i+1]]]
    l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[i+1]]]/M)
    w = 1/(1+l*GCovEst11[i,i])^2
    l = l[ind]
    z = z[ind]
    w = w[ind]
    z[which(z>zsquare_thresh)] = zsquare_thresh
    fit0=lm(z~l-1,weights=w)
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
  }
  
  lower_bounds <- c(Boundary$ecov.lower, Boundary$gcov.lower)
  upper_bounds <- c(Boundary$ecov.upper, Boundary$gcov.upper)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      ind = union(which(ZMatrix1[[col_names[i+1]]]!=0),which(ZMatrix1[[col_names[j+1]]]!=0))
      M1 = length(M1)
      z = ZMatrix1[[col_names[i+1]]] * ZMatrix1[[col_names[j+1]]]
      l = LDSC1$LDSC * sqrt(NMatrix1[[col_names[i+1]]]/M) * sqrt(NMatrix1[[col_names[j+1]]]/M)
      li = LDSC1$LDSC * NMatrix1[[col_names[i+1]]]/M
      lj = LDSC1$LDSC * NMatrix1[[col_names[j+1]]]/M
      w=1/((1+li*GCovEst[i,i])*(1+lj*GCovEst[j,j])+(l*GCovEst11[i,j]+ECovEst11[i,j])^2)
      z = z[ind]
      l = l[ind]
      li = li[ind]
      lj = lj[ind]
      w = w[ind]
      z[which(abs(z)>cov_thresh)] = sign(z[which(abs(z)>cov_thresh)]) * cov_thresh
      z1=z-ECovEst11[i,j]
      fit0=lm(z1~l-1,weights=w)
      X=cbind(1,l)
      a=c(sqrt(ECovEst[i,i]*ECovEst[j,j]),sqrt(GCovEst[i,i]*GCovEst[j,j]))
      result <- nloptr(x0 = c(ECovEst11[i,j],GCovEst11[i,j])*0.5,
                       eval_f = objective,
                       lb = lower_bounds*a,
                       ub = upper_bounds*a,
                       opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
      beta=result$solution
      vare=var(c(z-X%*%beta)*sqrt(w))
      H=solve(t(X)%*%(X*sqrt(w)))*vare
      GCovEst[i,j]=GCovEst[j,i]=beta[2]
      ECovEst[i,j]=ECovEst[j,i]=beta[1]
      GCovSE[i,j]=GCovSE[j,i]=sqrt(H[2,2])
      ECovSE[i,j]=ECovSE[j,i]=sqrt(H[1,1])
    }
  }
  
  t2=Sys.time()-t2
  print("Final Genetic Covariance Estimate")
  print(t2)
  
  GCovEst=positive_adj(GCovEst,min.eps=min.eps)
  ECovEst=positive_adj(ECovEst,min.eps=min.eps)
  row.names(GCovEst)=colnames(GCovEst)=row.names(ECovEst)=colnames(ECovEst)=NAM
  row.names(GCovSE)=colnames(GCovSE)=row.names(ECovSE)=colnames(ECovSE)=NAM
  row.names(GCovEst1)=colnames(GCovEst1)=row.names(ECovEst1)=colnames(ECovEst1)=NAM
  row.names(GCovSE1)=colnames(GCovSE1)=row.names(ECovSE1)=colnames(ECovSE1)=NAM
  Estimate.ini=list(GCovEst=GCovEst1,GCovSE=GCovSE1,ECovEst=ECovEst1,ECovSE=ECovEst1)
  Computing.time=list(stage1.time=t0,stage2.time=t1,stage3.time=t2)
  return(A=list(GCovEst=GCovEst,GCovSE=GCovSE,ECovEst=ECovEst,ECovSE=ECovEst,Estimate.ini=Estimate.ini,Computing.time=Computing.time))
}

