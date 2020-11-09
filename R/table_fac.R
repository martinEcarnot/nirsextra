table_fac = function (obs,pred,fac) {

cat("\nGlobal\n")
print(table(Classes_Obs=obs,Classes_pred=pred))
# ls=lapply(seg,FUN=function(x) {lapply(x,FUN=sort)}) # renge les indices de seg ds le mm ordre que la sortie de fitcv

nfac=ncol(fac)
for (i in 1:nfac) {
  fac1=pull(fac,i)
  for (j in 1:nlevels(fac1)){
    lev1=fac1==levels(fac1)[j]
    cat("\nFacteur ",colnames(fac)[i]," : ",levels(fac1)[j]," \n")
    print(table(Classes_Obs=obs[lev1],Classes_pred=pred[lev1]))
  }
}
}
