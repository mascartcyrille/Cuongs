#include "lassoshooting.hpp"
mat lasso_shooting(const mat &G, const umat & b, const mat & d, const vec & a0, int nitmax, double eps, uvec & nit)
{
  size_t M = b.n_cols;
  size_t Dim = b.n_rows;  
  mat a(Dim, M, fill::zeros);
  nit.set_size(M); 
  nit.zeros();
  vec w;
//  cout << "G = " << endl;
//  cout << G << endl;
//  cout << endl;
//  cout << "b = " << endl;
//  cout << b << endl;
//  cout << endl;
//  cout << "d = " << endl;
//  cout << d << endl;
//  cout << endl;
//  cout << "a0 = " << endl;
//  cout << a0 << endl;
//  cout << endl;
  for(size_t r=0; r<M; r++) {
    //    lasso_shooting1(G, b.col(r), d.col(r), a0, nitmax, eps, w, nit(r));
    lasso_shooting1(G, b.submat(0, r, Dim-1, r), d.submat(0, r, Dim-1, r), a0, nitmax, eps, w, &nit(r));
    a.submat(0, r, Dim-1, r) = w;
  }
  return a;
}

double J(const mat &G, const uvec & b, const vec & d, const vec & a)
{
  return 1./2*dot(a, G*a) - dot(b, a) + dot(d, arma::abs(a));
}

void lasso_shooting1(const mat &G, const uvec & b, const mat & d, const vec & a0, int nitmax, double eps, vec & a, long long unsigned int*  nit)
//void lasso_shooting1(const mat &G, const uvec & b, const mat & d, const vec & a0, int nitmax, double eps, vec & a, size_t nit)
{
  size_t Dim = b.size();
  vec aold;
  double S;
  a = a0;
  aold = a;
  //  nit = 0;
  *nit = 0;
  //  while((( std::abs(J(G, b, d, a) - J(G, b, d, aold))> eps) & (nit<=nitmax))| (nit==0) )
  while((( std::abs(J(G, b, d, a) - J(G, b, d, aold))> eps) & ((*nit)<=nitmax))| ((*nit)==0) )
    {
      aold = a;
      for(size_t i=0; i<Dim; i++)
	{
	  S = dot(G.row(i),a) - G(i,i)*a(i) - b(i);
	  if(G(i,i)<=0.)
	    cout << "Unexpected value on the diagonal of G'" << endl;
	    //	    warning('Unexpected value on the diagonal of G');
	  if((-S-d(i))>eps) 
	    a(i) = (-S - d(i))/G(i,i);
	  else if ((-S+d(i) )< -eps)
	    a(i) = (-S + d(i))/G(i,i);
	  else
	    a(i) = 0;
	}
      //      nit++;
      (*nit)++;
    }
  //cout << "nit = " << nit << endl;
  //cout << "nit = " << nit << endl;
}
