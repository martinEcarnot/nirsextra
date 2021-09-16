table_fac = function (obs,pred,fac) {

# library(dplyr)

# Fonction pour afficher le % de bien classe dans une matrice de mm dim que la table de contingence
make_pbc= function(tabl,nlev) {
pbc=100*sum(diag(tabl))/sum(tabl)
t=array(c(c(sprintf("  %.1f%%  ",pbc),rep("",nlev-1)),rep("",nlev)),dim = c(nlev,nlev))
return(t)
}

obs=as.factor(obs)
pred=as.factor(pred)

nlev=nlevels(as.factor(obs))
cat("\nGlobal\n")
t=table(obs,pred)
disp=cbind(t,make_pbc(t,nlev))
rownames(disp)=c("Classes_Obs",rep("",nlev-1))
colnames(disp)=c("Classes_pred","","","")
print(disp,quote=F)

# ls=lapply(seg,FUN=function(x) {lapply(x,FUN=sort)}) # renge les indices de seg ds le mm ordre que la sortie de fitcv

nfac=ncol(fac)
for (i in 1:nfac) {
  fac1=pull(fac,i)
  disp=matrix(nrow=nlev,ncol=0)
  cname=character()
  for (j in 1:nlevels(fac1)){

    lev1=fac1==levels(fac1)[j]
    t=table(Classes_Obs=obs[lev1],Classes_pred=pred[lev1])
    disp=cbind(disp,t,make_pbc(t,nlev))
    cname=c(cname,"Pred","",levels(fac1)[j],"")  # c(cname,"Classes_pred","",levels(fac1)[j],"")
  }
  # browser()
  rownames(disp)=c("Obs",rep("",nlev-1))  # c("Classes_Obs",rep("",nlev-1))
  colnames(disp)=cname #rep(c("Classes_pred","","",""),nlevels(fac1))
  cat("\n\nFacteur: ",colnames(fac)[i]," \n\n")
  print(disp,quote=F)
}


}



