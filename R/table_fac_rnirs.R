table_fac_rnirs = function (fm,dat_fac,numfac,segm,comp,rep=1:max(fm$y$rep)) {

  fm = onefitcv(fm,comp,rep)

  table_fac(fm$y$y1,fm$fit$y1,dat_fac[unlist(segm),numfac])
}



