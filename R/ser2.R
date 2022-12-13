SER2 <- function(yp,yo)
{
  R2=cor(yp,yo)^2
  if (is.null(ncol(yp))) {
    SE=sqrt(mean((yp-yo)^2))
  } else {
    SE=sqrt(colMeans((replicate(ncol(yp),yo)-yp)^2))
  }
  ser2=data.frame(SE,R2)
  return(ser2)
}
