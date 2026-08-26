#' adj_NaturaSpec
#'
#'  # Elimine les sauts dus au passage d'un dtecteur  a un autre.
#'  On extrapole l'ajustement linaire fait sur les [ws] points avant le saut, sur les 14 points suivants
#'  (interpolated after gap, I guess from SpectralEvolution)
#'  La difference entre le point à gap + 15 (entre vrai spectre et fit) est soustraite a tous les points suivants
#'  [632,?]
#'
#' @export
#'

adj_NaturaSpec <- function(X) {

  ws=5  # Number of wavelengths used for linear adjustment before gap
  interpw=14  # Number of wavelengths interpolated after gap (I guess from SpectralEvolution)

  for (iadj in c(632)) {

    x= (iadj - ws + 1):iadj
    Y=Xi[, x]
    my=t(colMeans(t(Y)))

    sx=var(x)
    mx=mean(x)
    b=cov(x,t(Y))/sx # Ajustement linaire.
    b0=my-b*mx  # Ajustement linaire.

    # 14 points after 632 seems to be interpolated: replaced by linear adjustment from 5 last points
    Xi[,(iadj+1):(iadj+interpw)] = drop(replicate(interpw,b0)) + t(b)%*%((iadj+1):(iadj+interpw))

    # Compute difference between adjusted and original spectra after gap, and apply it to the rest of the spectrum
    dif=Xi[,iadj+interpw+1]-(b0+b*(iadj+interpw+1))
    loremp=(iadj+interpw+1):dim(Xi)[2]
    Xi[,loremp]=Xi[,loremp] - drop(replicate(length(loremp),dif))
  }
  return(Xi)

}
