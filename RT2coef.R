RT2coef <-function(fm)  {

# From fm model (rnirs/rchemo) computes Beta coef


b=fm$R[,1:i] %*% ginv(t(fm$Tr[,1:i]) %*% fm$Tr[,1:i]) %*% t(fm$Tr[,1:i]) %*% yarn$density
}
