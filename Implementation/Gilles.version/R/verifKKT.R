verifKKT <- function(x, G, b, d, eps=1e-7)
{
  v1 = G%*%x - b + d*sign(x)
  Ip = which(abs(x)>eps)
  Iz = which(abs(x)<=eps)
  if(is.matrix(x)) {
    aux = t(x)%*%G%*%x-t(x)%*%b+ d%*%abs(x) }
  else {
    aux = x%*%G%*%x-x%*%b+ t(abs(x))%*%d }
  
  print(paste('Checking x^tGx - b^tx + d^t|x| should be zero = ', aux))
  return(list(all(abs(v1[Ip]<=eps)) * all((v1[Iz]>= -d [Iz]) & (v1[Iz]<= d [Iz])), aux))
}
