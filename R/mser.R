#' mser
#'
#' Computes statistics to compare prediction and reference values
#'
#' @import dplyr
#' @export

mser = function(fm)  {

y=fm$y

  res <- y %>%
    group_by(nlv) %>%
    summarise(
      npred = n(),
      sep = sep(yp,yref),
      cor2 = cor2(yp,yref),
      RDP = rpd(yp,yref),
      .groups = "drop"
    ) %>%
    arrange(nlv)

}
