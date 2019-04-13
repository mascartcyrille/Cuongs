#include "lassoshooting.hpp"

mat set_d(int M, int K, double gamma, const umat & mu2, const umat & muA)
{
  int Dim  = mu2.n_elem;
  double mylog = std::log((1+M*K)*M);
  mat P = sqrt(2*gamma*mylog*conv_to<arma::mat>::from(mu2)) + gamma/3.*mylog*conv_to<arma::mat>::from(repmat(muA, 1, M));
  return P;
}

int main(void){
  mat G = {{10, 5, 3, 2, 1},{5, 12, 4, -1, -1}, {3, 4, 15, 6, -1}, {2, -1, 6, 17, 3}, {1, -1, -1, 3, 20}};
  double gamma = 3;
  int M, K;
  M = K = 2;
  umat mu1 = {{5, 6}, {19, 2}, {7, 25}, {4, 2}, {10, 6}};
  umat mu2 = {{6, 17}, {12, 2}, {4, 2}, {4, 8}, {6, 9}};
  uvec muA = {2, 1, 1, 3, 7};
  int n;
  mat a;
  vec a0(5, fill::ones);
  a0 = -17*a0;
  mat d;
  d = set_d(M, K, gamma, mu2, muA);
  uvec nit;
  double eps = 1e-12;
  int nitmax = 10000;
  a = lasso_shooting(G, mu1, d, a0, nitmax, eps, nit);
  //  cout.precision(11);
  //  cout.setf(ios::fixed);
  cout<< "Our Lasso solution is "<< endl << a << endl;
  cout << "Number of iterations " << endl << nit << endl;
  cout << "R lassoshooting package gives the following solution: [0, 0.3184682, 0, 0, 0], [0, 0, 1.162726, 0, 0] (see ../R/test_R_lasso.R file)" << endl;
}
