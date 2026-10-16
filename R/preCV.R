preCV = function(x,y,fun="plskern",list_pre,ncomp,titl, plotLV=TRUE, plotYY=TRUE, verb = FALSE)  {

  # From a list of pretraitments, make cross-validation of data set with different pretr. and plot pred vs obs for best LV

  r2_tt=matrix(nrow=ncomp+1,ncol=length(list_pre))  # matrix(nrow=ncomp+1,ncol=nseq-6)  #
  fmtt=list()
  pftot=NULL
  K=5
  # nseq=20   # For testing CV on slices of spectra
  # seqlo=seq(1,2150,nseq)
  # nseq=length(seqlo)
  r2_tt=matrix(nrow=ncomp,ncol=length(list_pre))  # matrix(nrow=ncomp+1,ncol=nseq-6)  #
  fmtt=list()
  fmttn=list()
  for (j in 1:length(list_pre)) {  # 1:
    # print(paste0("Prétraitement ",j, "/",length(list_pre),"\n"))
    if (verb) {cat(paste0("Prétraitement ", j, "/", length(list_pre), "\n"))}
    # p=rbind(list('adj',''),list('red',c(seqlo[j],2151-seqlo[j+6],1)),list('sder',c(1,3,15)))
    # print(seqlo[j])
    # print(2151-seqlo[j+4])
    # dat$xp=pre(dat$x,p)
    xp=pre(x,list_pre[[j]])

    # generate sgm list for Leave-one-out
    # segm <- list(rep1 = as.list(1:nrow(xp)))
    segm <- segmkf(n = nrow(xp), K = K)
    if (fun=="lwplsr") {
      pars=mpars(nlvdis = 5, diss = "mahal", h = c(1, Inf), k = c(10, 20))
      fmc = gcvlv(xp, y,segm,score = r2, fun = get(fun), pars=pars, nlv = 1:ncomp, verb = F)  # !!! pas cor2 avec LOO
    }
    else {
      fmc = gcvlv(xp, y,segm,score = r2, fun = get(fun), nlv = 1:ncomp, verb = F)  # !!! pas cor2 avec LOO
    }
    fmtt=append(fmtt, list(fmc)) # fmttn=append(fmttn, list(fm))
    r2_tt[,j] = mser(fmc)$cor2  # r2_tt[,j] = mse(fm, ~ ncomp)$cor2
    pf=t(list_pre[[j]])  # # list(seqlo[j], 2151-seqlo[j+4])
    dim(pf)=c(1,2*ncol(pf))
    pf=paste0(pf, collapse = "_")
    pftot=c(pftot,pf)
  }
  best_pre_lo=which(r2_tt == max(r2_tt[2:ncomp,]), arr.ind = TRUE)[1,]

  if (plotLV) {
    matplot(r2_tt, type = 'l', lty = 1, col = 1:ncol(r2_tt), ylab="R2_Cross-Validation", xlab="Number of Latent Variables")
    title(titl)
    legend("bottomright", legend = pftot, col = 1:ncol(r2_tt), lty = 1, cex = 0.8,bg = "white")
  }
  fm=fmtt[[best_pre_lo[2]]]  # fmn=fmttn[[best_pre_lo[2]]]
  txtpre= pftot[best_pre_lo[2]]
  cat(
    "Pré :", txtpre, "\n\n",
    "R2     :", round(mser(fm)$cor2, 2), "\n\n",
    "SEP    :", round(mser(fm)$sep, 2), "\n\n\n"
  )

  if (plotYY) {
    fm1=fm$y[fm$y$nlv==best_pre_lo[1],]
    plot(fm1$yp,fm1$yref, xlab="NIRS prediction (CV)", ylab="Reference value")
    abline(a = 0, b = 1, col = "red", lty = 2)  # a=0, b=1 pour y=x; couleur rouge et ligne en pointillés
    fit <- lm(fm1$yp ~ fm1$yref)
    abline(fit, col = "blue")
    summary_fit <- summary(fit)
    r_squared <- summary_fit$r.squared
    secv=sep(fm1$yp,fm1$yref)
    legend("topleft", legend = c("y = x", bquote(Cross_Validation: ~ R^2 == .(round(r_squared, 2))),bquote(SECV : .(round(secv, 2))),bquote(ncomp : .(best_pre_lo[[1]])),bquote(n_smpl : .(length(fm1$yref)))), col = c("red", "blue", "white", "white", "white"),lty = c(2, 1),bty = "n")
    title(titl)
      }
}


# disp_pre=function(x,pre1) {
#   txt=NULL
#   for (i in 1:nrow(pre1)) {
#     if (pre1[i,1] == "adj") {}
#     if (pre1[i,1] == "snv") {txt=paste0(txt,"SNV")}
#     if (pre1[i,1] == "red") {
#       txt=paste0(txt,"RANGE:",min(as.numeric(attributes(x)[2]$dimnames[[2]])), "-",max(as.numeric(attributes(x)[2]$dimnames[[2]])),"nm")
#     }
#     if (pre1[i,1] == "sder") {
#       if (pre1[i,2][[1]][1] == 1) {txtorder= "1ere"}
#       if (pre1[i,2][[1]][1] == 2) {txtorder= "2nde"}
#       txt=paste0(txt,"DERIV.",txtorder, " (","poly=",pre1[i,2][[1]][2],";","window=",pre1[i,2][[1]][2],")")
#     }
#     if (length(txt)>0 & i<nrow(pre1)) {txt=paste0(txt," + ")}
#   }
#   return(txt)
# }
