block.generate=function(M,nblock){
start=seq(1,M-floor(M/nblock),floor(M/nblock))
end=start+floor(M/nblock)-1
return(list(start=start,end=end))
}
