#include "misc.hpp"
#include "DataSpike.hpp"

//MatrixXd computeGk1k2_V2(int M, int K, double Tmin, double Tmax, VectorXd& T, VectorXi & neur, double delta, VectorXi&low)
mat computeGk1k2_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec &low)
{
  const vec & T = DS._T;
  const uvec & neur = DS._neur;
  mat G((1+M*K), (1+M*K), fill::zeros);
   G(0,0) = Tmax - Tmin;
  double A = K*delta;
  int depart = get_low_index(Tmin-A, T);
  int beta = get_low_index(Tmax, T) - 1; // get_up_index(Tmax, T);
  double ti, tj, x1, x2, dx;
  int l1, l2;
  umat Mk;
  size_t k1, k2;
  for(int i=depart;i<beta+1;i++)
    {
      ti = T(i); l1 = neur(i);
      // First row and column
      for(int k=1;k<=K;k++)
	{
	  x1 = min(Tmax, ti+k*delta);
	  x2 = max(Tmin, ti+(k-1)*delta);
	  dx = x1 - x2;
	  if(dx >0)
	    {
	      G(0,(l1-1)*K+k) += dx;
	      G((l1-1)*K+k,0) += dx;
	    }
	}
      for(int j=low(i);j<i;j++)
	{
	  tj = T(j); l2 = neur(j);
	  Mk = get_k1k2(K, ti - tj, delta); 
	  //	  cout << "Mk.n_rows " << Mk.n_rows << endl;
	  for(size_t iel=0; iel<Mk.n_rows; iel++)
	    {
	      k1 = Mk(iel,0); k2 = Mk(iel,1);
	      x1 = min(Tmax, min(ti+k1*delta, tj+k2*delta));
	      x2 = max(Tmin, max(ti+(k1-1)*delta, tj+(k2-1)*delta));
	      dx = x1 - x2;
	      if(dx >0)
		{
		  G((l1-1)*K+k1,(l2-1)*K+k2) += dx;
		  G((l2-1)*K+k2,(l1-1)*K+k1) += dx;
		}
	    }
	}
      
      for(int k1=1;k1<=K;k1++)
	{
	  x1 = min(Tmax, ti+k1*delta);
	  x2 = max(Tmin, ti+(k1-1)*delta);
	  dx = x1 - x2;
	  if(dx >0)
	    G((l1-1)*K+k1,(l1-1)*K+k1) += dx;
	  for(int k2=k1+1;k2<=K;k2++)
	    {
	      x2 = max(Tmin, ti+(k2-1)*delta);
	      dx = x1 - x2;
	      if(dx >0)
		{
		  G((l1-1)*K+k1,(l1-1)*K+k2) += dx;
		  G((l1-1)*K+k2,(l1-1)*K+k1) += dx;
		}
	    }
	  
	}
    }
  return G;
}

mat computeGk1k2_V2(int M, int K, double Tmin, double Tmax, const mat & DN, double delta, uvec &low)
{
  DataSpike DS(DN);
  mat G = computeGk1k2_V2(M, K, Tmin, Tmax, DS, delta, low); 
  return G;
}
