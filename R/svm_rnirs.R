#' svm_rnirs
#'
#' svm_rnirs
#'
#' @export
#'

svm_rnirs <- function (Xr,Yr,Xu,Yu,cost=NULL, epsilon=NULL, gamma=NULL)  {


  # Copié de lmr.R
  # Xr <- as.matrix(Xr, prefix.colnam = "x")
  n <- nrow(Xr)
  p <- ncol(Xr)
  # Xu <- as.matrix(Xu, prefix.colnam = "x")
  m <- nrow(Xu)
  rownam.Xu <- row.names(Xu)
  # Yr <- as.matrix(Yr, row = FALSE, prefix.colnam = "y")
  q <- ncol(Yr)
  colnam.Yu <- colnames(Yr)
  if (is.null(Yu))
    Yu <- matrix(nrow = m, ncol = q)
  else Yu <- as.matrix(Yu, row = FALSE, prefix.colnam = "y")

  xcal=sp2df(Xr,Yr)   # ajout
  xval=sp2df(Xu,Yu)   # ajout

  if (is.null(cost))
    cost <- 1
  if (is.null(epsilon))
    epsilon <- 0.1
  if (is.null(gamma))
    gamma <- 1 / ncol(Xr)

  model <- svm(y ~ x , xcal,cost=cost , epsilon=epsilon, gamma=gamma)   # ajout


  y <- Yu
  fit <- predict(model, xval)   # ajout
  r <- y - fit
  dat <- data.frame(rownum = 1:m, rownam = rownam.Xu)
  y <- cbind(dat, y)
  fit <- cbind(dat, fit)
  r <- cbind(dat, r)
  zq <- ncol(y)
  u <- (zq - q + 1):zq
  names(r)[u] <- names(fit)[u] <- names(y)[u] <- colnam.Yu
  list(y = y, fit = fit, r = r)



}
