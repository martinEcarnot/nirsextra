#' naturaspec2df
#'
#' Read NaturaSpec files and set it to a dataframe
#'
#' @export
#'
naturaspec2df <- function(d) {
  sp=read_spectra(d)
  x=sp$value
  colnames(x)=sp$bands
  row.names(x)=substr(sp$names,1,nchar(sp$names)-10)
  df=sp2df(x)
}
