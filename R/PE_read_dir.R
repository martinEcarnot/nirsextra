PE_read_dir <- function(d,nlo)
{
  # Reading PE files on a directory
  #

  if (missing(nlo)) {nlo=3001}
  l=Sys.glob(file.path(d,"*.sp"))
  x=matrix(, nrow = length(l), ncol = nlo)

  for (i in 1:length(l)) {
    print(l[i])
    sp=readSP.pe(l[i])
    x[i,]=sp$yData
  }
  row.names(x)=gsub(".sp","",basename(l))
  colnames(x)=sp$xData
  return(x)
}


