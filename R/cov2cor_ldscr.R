cov2cor_ldscr=function(X){
x=diag(X)
x[x<(1e-10)]=0
x[x>0]=1/x[x>0]
x=sqrt(x)
B=t(t(X)*x)*x
return(B)
}
