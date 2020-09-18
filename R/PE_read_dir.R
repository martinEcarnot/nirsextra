PE_read_dir <- function(d)
{
  # Reading PE files on a directory
  #

  l=Sys.glob(file.path(d,"*.sp"))
  x=matrix(, nrow = length(l), ncol = 3001)

  for (i in 1:length(l)) {
    sp=readSP.pe(l[i])
    x[i,]=sp$yData
  }
  row.names(x)=gsub(".sp","",basename(l))
  colnames(x)=sp$xData
  return(x)
}


