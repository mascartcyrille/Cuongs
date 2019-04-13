source('verifKKT.R')
library(lassoshooting)

                                        # definition of matrix D and vector x from G, b and d (D and x will be used in R lassoshooting library)
                                        # returns a list containing D and x
#
lasso_convert<- function(G, b, d)
{
  eps = 1e-8
  U = chol(G)
  x = solve(t(U), b)
  stopifnot(all(abs(d)>eps))
  D = U %*% diag(1./d)
  return(list(D, x))
}

                                        # definition of d from mu2, muA, gamma, M, K
#
set_d <- function(M, K, gamma, mu2, muA)
{
  return(sqrt(2*gamma*log((1+M*K)*M)*mu2)+gamma/3.*log((1+M*K)*M)*matrix(rep(muA, M), nrow= length(muA)))
}

# example
G = matrix(c(10, 5, 3, 2, 1, 5, 12, 4, -1, -1, 3, 4, 15, 6, -1, 2, -1, 6, 17, 3, 1, -1, -1, 3, 20), nrow=5)
M = K = 2
gamma = 3
mu1 = matrix(c(5, 19, 7, 4, 10, 6, 2, 25, 2, 6), ncol =2)
mu2 = matrix(c(6, 12, 4, 4, 6, 17, 2, 2, 8, 9), ncol=2)
muA = c(2, 1, 1, 3, 7)
d = set_d(M, K, gamma, mu2, muA)
for(r in seq_len(M)) {
  L = lasso_convert(G, mu1[,r], d[,r])
  x = lassoshooting(L[[1]], L[[2]], 1.)
  a = diag(1./d[,r])%*% x$coefficients
  print(paste('Lasso solution for neuron', r))
  print(a)
  print(verifKKT(a, G, mu1[,r], d[,r])) # verification of KKT conditions (to check if the solution is correct from a mathematical point of view)
}
