#' orthEf
#'
#' Je sais meme plus ce que ça fait...
#'
#' @export

orthEf <- function (datref,datappli,spname,coleffet,ncomp) {

# ieffet=which(colnames(datref) %in% effet)
S=datref[,which(colnames(datref) %in% spname)]

  Ds=tab.disjonctif(datref[,coleffet])
  # browser()
  BS =  Ds %*% pinv(t(Ds) %*% Ds)$Xplus %*% t(Ds) %*% S
  WS = S - BS;

# ACP
  rpca=pcasvd(WS,ncomp= ncomp)  # rpca=pca(WS,ncomp= qr(WS)$rank)
loadings_error=rpca$P
kW= loadings_error[,1:ncomp]

datappli$xort=datappli[,spname] - datappli[,spname] %*% kW %*% t(kW)

  return(datappli)
}
