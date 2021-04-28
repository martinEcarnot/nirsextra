asd_read_dir <- function(d,times=F, comments=F)
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
  x[i,]=sp$spectrum/sp$reference
  # x[i,]=sp$spectrum
  if (times) {
  t1=paste(sp$info$spectrumheaderwhensec,sp$info$spectrumheaderwhenmin,sp$info$spectrumheaderwhenhour,sp$info$spectrumheaderwhenmday,sp$info$spectrumheaderwhenmon,sp$info$spectrumheaderwhenyear)
  tim=rbind(tim,as.POSIXct(t1, format = "%S %M %H %d %m %y",origin = "1960-01-01"))
  }
  if (comments) {com=rbind(com,sp$info$spectrumheadercomments)}
}
if (times) {attr(x,"time")=tim[-1]}
  row.names(x)=gsub(".asd","",basename(l))
  colnames(x)=seq(from=sp$info$spectrumheaderch1_wavel,length.out=sp$info$spectrumheaderwchannels, by=sp$info$spectrumheaderwavel_step)
if (comments) {
  x=sp2df(x)
  x$com=com
}


return(x)
}


