#' Single-Variate Linkage Disequilibrium Score Regression (LDSC)
#'
#' The `ldsc.uni` function performs single-variate Linkage Disequilibrium Score Regression (LDSC) analysis. It is designed to estimate heritability from GWAS summary statistics, accounting for linkage disequilibrium (LD) between SNPs. The function harmonizes GWAS data with LD scores and applies non-linear optimization to estimate heritability.
#'
#' @param gwas A data.frame containing GWAS summary statistics for a single trait. The data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param LDSC A data.frame containing LD Score Regression (LDSC) estimates. It should include LDSC scores and other necessary metrics for the analysis.
#' @param Boundary An optional list of parameter estimate boundaries. If not provided, default values are used. It typically includes intercept.lower (the lower boundary of the intercept estimate), intercept.upper (the upper boundary of the intercept estimate), and h2.upper (the upper boundary of the heritability estimate).
#' @param zsquare_thresh A threshold for the squared Z-scores in heritability estimation to control for extreme values.
#' @param min.eps A small positive value to ensure numerical stability in calculations.
#' @param nblock The number of blocks for bootstrap-based standard error estimation.
#' @param sampling.time The number of block bootstrap.
#' @param sampling.ratio The sub-sampling ratio in each bootstrap.
#' @return A data.frame containing heritability estimates and their standard errors, along with the intercept and its standard error.
#'
#' @examples
#' # Example usage
#' data(hapmap3)
#' data(EURLDSC)
#' ref_panel <- hapmap3
#' LDSC <- EURLDSC
#' results <- ldsc.univ(gwas = GWAS_data, LDSC = LDSC)
#'
#' @details The `ldsc.univ` function is particularly useful for single-trait heritability analysis using LDSC. It requires well-formatted GWAS summary statistics and an understanding of LDSC methodology. The function includes an initial estimation phase followed by a reweighting step for improved efficiency and accuracy.
#'
#' @importFrom stats lm
#' @importFrom data.table setDT setkey copy setnames
#' @importFrom nloptr nloptr
#'
#' @export
#'
ldsc.univ=function(gwas,LDSC,Boundary=F,zsquare_thresh=500,nblock=200,sampling.time=200,sampling.ratio=0.5){

  ############################# Basic Information ###############################
  t0 = Sys.time()
  gwas=merge(gwas,LDSC,by="SNP")
  M=nrow(gwas)
  if(Boundary[1]==F){
    Boundary=list(intercept.lower=0.95,intercept.upper=1.05,h2.upper=0.95)
  }
  t0 = Sys.time() - t0
  print("Processing data")
  print(t0)

  ############################# initial estimator ####################################
  t1=Sys.time()
  objective <- function(beta) {
    sum((z - X %*% beta)^2)
  }
  lower_bounds <- c(Boundary$intercept.lower, 0)
  upper_bounds <- c(Boundary$intercept.upper, Boundary$h2.upper)

  z = gwas$Zscore^2
  l = gwas$LDSC*gwas$N/M
  z[which(z>zsquare_thresh)] = zsquare_thresh
  X=cbind(1,l)
  result <- nloptr(x0 = c(1,0.1),
                   eval_f = objective,
                   lb = lower_bounds,
                   ub = upper_bounds,
                   opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
  beta=result$solution
  h2.ini=beta[2]
  intercept.ini=beta[1]

  t1=Sys.time()-t1
  print("Initial Genetic Covariance Estimate")
  print(t1)

  ############################ reweight for efficiency ##################################
  t2=Sys.time()
  objective <- function(beta) {
    sum((z - X %*% beta)^2*w)
  }
  z = gwas$Zscore^2
  l = gwas$LDSC*gwas$N/M
  z[which(z>zsquare_thresh)] = zsquare_thresh
  w = 1/(1+l*h2.ini)^2
  X=cbind(1,l)
  result <- nloptr(x0 = c(intercept.ini,h2.ini),
                     eval_f = objective,
                     lb = lower_bounds,
                     ub = upper_bounds,
                     opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
  beta=result$solution
  h2=beta[2]
  intercept=beta[1]
  t2=Sys.time()-t2
  print("Heritability Estimation")
  print(t2)

  t3=Sys.time()
  objective <- function(beta) {
    sum((z1 - X1 %*% beta)^2*w1)
  }
  w = 1/(1+l*h2)^2
  blockind=block.generate(M,nblock)
  h2.vec=intercept.vec=c(1:sampling.time)
  for(i in 1:sampling.time){
  indblock=sample(nblock,nblock*sampling.ratio)
  ind=NULL
  for(j in 1:length(indblock)){
  ind=c(ind,c(blockind$start[indblock[j]]:blockind$end[indblock[j]]))
  }
  z1=z[ind]
  l1=l[ind]
  w1=w[ind]
  X1=cbind(1,l1)
  result <- nloptr(x0 = c(intercept,h2),
                   eval_f = objective,
                   lb = lower_bounds,
                   ub = upper_bounds,
                   opts = list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=100))
  beta=result$solution
  h2.vec[i]=beta[2]
  intercept.vec[i]=beta[1]
  }
  t3=Sys.time()-t3
  print("Standard Error Estimation")
  print(t3)
  intercept.se=sqrt(mean((intercept.vec-intercept)^2))
  h2.se=sqrt(mean((h2.vec-h2)^2))
  return(A=data.frame(intercept=intercept,intercept.se=intercept.se,h2=h2,h2.se=h2.se))
}

