onefitcv = function (fmncmp,ncomp,rep=1:max(fm$y$rep)) {

  # Pour extraire une compsante de fm genere par fitcv
    fm=lapply(fmncmp[1:3],function (x) {x[x$ncomp==ncomp & x$rep %in% rep,]})
}
