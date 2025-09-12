#' asd_read_dir
#'
#' Read all .asd file from repertoty
#'
#' @export
#'

asd_read_dir <- function(d,times=F, comments=F, type=F)
{
  # Reading ASD files on a directory
  #

l=Sys.glob(file.path(d,"*.asd"))
x=matrix(, nrow = length(l), ncol = 2151)
tim=""
com=NULL

for (i in 1:length(l)) {
  # x[i,]=get_spectra(l[i], type = "reflectance")
  sp=asd_read(l[i])
  if (type=="l") {
    x[i,]=sp$spectrum
  } else if(type=="ref") {
    x[i,]=sp$reference
  } else {
    x[i,]=sp$spectrum/sp$reference
  }
  # x[i,]=sp$spectrum
  if (times) {
  t1=paste(sp$info$spectrumheaderwhensec,sp$info$spectrumheaderwhenmin,sp$info$spectrumheaderwhenhour,sp$info$spectrumheaderwhenmday,sp$info$spectrumheaderwhenmon,sp$info$spectrumheaderwhenyear)
  # tim=rbind(tim,as.POSIXct(t1, format = "%S %M %H %d %m %Y"))
  tim=rbind(tim,t1)
  }
  if (comments) {com=rbind(com,sp$info$spectrumheadercomments)}
}
# if (times) {attr(x,"time")=tim[-1]}
if (times) {attr(x,"time")=as.POSIXct(tim[-1], format = "%S %M %H %d %m %Y")}
  row.names(x)=gsub(".asd","",basename(l))
  colnames(x)=seq(from=sp$info$spectrumheaderch1_wavel,length.out=sp$info$spectrumheaderwchannels, by=sp$info$spectrumheaderwavel_step)
if (comments) {
  x=sp2df(x)
  x$com=com
}


return(x)
}


