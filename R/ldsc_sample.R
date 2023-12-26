ldsc_sample=function(GWAS_List,LDSC){
############################# Basic Information ###############################
NAM=names(GWAS_List)
p=length(GWAS_List)
ZMatrix=data.frame(SNP=GWAS_List[[1]]$SNP)
for(i in 1:p){
A=GWAS_List[[i]]
ZMatrix=cbind(ZMatrix,A$Zscore)
}
names(ZMatrix)=c("SNP",NAM)
NMatrix=data.frame(SNP=GWAS_List[[1]]$SNP)
for(i in 1:p){
A=GWAS_List[[i]]
NMatrix=cbind(NMatrix,A$N)
}
names(NMatrix)=c("SNP",NAM)
row.names(NMatrix)=row.names(ZMatrix)=ZMatrix$SNP
SNPInfo=GWAS_List[[1]][,c("SNP","CHR")]

snplist=intersect(SNPInfo$SNP,LDSC$SNP)
ZMatrix1=ZMatrix[snplist,]
NMatrix1=NMatrix[snplist,]
LDSC1=LDSC[snplist,]
M=length(snplist)
############################# initial estimator ####################################
t1=Sys.time()
GCovEst1=GCovSE1=ECovEst1=ECovSE1=diag(p)*0
for(i in 1:p){
for(j in 1:p){
z=ZMatrix1[,i+1]*ZMatrix1[,j+1]
l=LDSC1$LDSC.R*sqrt(NMatrix1[,i+1]/M)*sqrt(NMatrix1[,i+1]/M)
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
for(i in 1:p){
z=ZMatrix1[,i+1]*ZMatrix1[,i+1]
l=LDSC1$LDSC.R*NMatrix1[,i+1]/M
w=(1+l*GCovEst[i,i])^2
fit0=lm(z~l,weights=1/w)
summary0=summary(fit0)
GCovEst[i,i]=summary0$coefficients[2,1]
ECovEst[i,i]=summary0$coefficients[1,1]
GCovSE[i,i]=summary0$coefficients[2,2]
ECovSE[i,i]=summary0$coefficients[1,2]
}

for(i in 1:(p-1)){
for(j in (i+1):p){
z=ZMatrix1[,i+1]*ZMatrix1[,j+1]
l=LDSC1$LDSC.R*sqrt(NMatrix1[,i+1]/M)*sqrt(NMatrix1[,i+1]/M)
li=LDSC1$LDSC.R*NMatrix1[,i+1]/M
lj=LDSC1$LDSC.R*NMatrix1[,j+1]/M
w=(1+li*GCovEst[i,i])*(1+lj*GCovEst[j,j])+(l*GCovEst1[i,j]+ECovEst1[i,j])^2
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
