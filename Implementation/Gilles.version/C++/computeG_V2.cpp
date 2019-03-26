#include "misc.hpp"
#include "DataSpike.hpp"

///
///\file computeG_V2.cpp
/// \brief Computation of the matrix G in the multivariate case
/// Algo proposed by PRB
///

//MatrixXd computeG_V2(int M, int K, double Tmin, double Tmax, VectorXd& T, VectorXi & neur, double delta, VectorXi&low)
mat computeG_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec &low)
{
///
/// \fn mat computeG_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec &low)
/// \brief Returns the matrix G of size (1+M*K)-by-(1+M*K)
/// \param M number of neurons
/// \param  K number of windows
/// \param  Tmin beginning time
/// \param  Tmax end time
/// \param  DS DataSpike instanceor like a matrix of size 2-by-Ntot, Ntot is the total number of spikes over all the neurons, the first row contains sorted spike times, the second one neuron numbers
/// \param  delta  K*delta is equal to the scope
/// \param  low array of size Ntot(= number of columns of DS) returned by computeb_V2 function, contains the first index where  
/// \return matrix G of size (1+M*K)-by-(1+M*K)
///
  const vec & T = DS._T;
  const uvec & neur = DS._neur;
  mat G((1+M*K), (1+M*K), fill::zeros);
  G(0,0) = Tmax - Tmin;
  double A = K*delta;
  // depart = integer part of Tmin-A
  // beta = integer part of Tmax
  int depart = get_low_index(Tmin-A, T);
  int beta = get_low_index(Tmax, T) - 1; 
  double ti, tj, x1, x2, dx;
  int l1, l2;
  for(int i=depart;i<beta+1;i++)
    {
      ti = T(i); l1 = neur(i);
      // First row and first column
      for(int k=1;k<K+1;k++)
	{
	  x1 = min(Tmax, ti+k*delta);
	  x2 = max(Tmin, ti+(k-1)*delta);
	  dx = x1 - x2;
	  if(dx >0)
	    {
	      G(0, (l1-1)*K+k) += dx;
	      G((l1-1)*K+k, 0) += dx;
	    }
	}
      for(int j=low(i);j<i;j++)
	{
	  tj = T(j); l2 = neur(j); 
	  for(int k1=1;k1<K+1;k1++)
	    {
	      for(int k2=1;k2<K+1;k2++)
		{
		  x1 = min(Tmax, min(ti+k1*delta, tj+k2*delta));
		  x2 = max(Tmin, max(ti+(k1-1)*delta, tj+(k2-1)*delta));
		  dx = x1 - x2;
		  if(dx >0)
		    {
		      G((l1-1)*K+k1, (l2-1)*K+k2) += dx;
		      G((l2-1)*K+k2, (l1-1)*K+k1) += dx;
		    }
		}
	    }
	}

      // diagonal part
      for(int k1=1;k1<K+1;k1++)
	{
	  x1 = min(Tmax, ti+k1*delta);
	  x2 = max(Tmin, ti+(k1-1)*delta);
	  dx = x1 - x2;
	  if(dx >0)
	    G((l1-1)*K+k1, (l1-1)*K+k1) += dx;
	  for(int k2=k1+1;k2<K+1;k2++)
	    {
	      x2 = max(Tmin, ti+(k2-1)*delta);
	      dx = x1 - x2;
	      if(dx >0)
		{
		  G((l1-1)*K+k1, (l1-1)*K+k2) += dx;
		  G((l1-1)*K+k2, (l1-1)*K+k1) += dx;
		}
	    }
	}
    }
  return G;
}

mat computeG_V2(int M, int K, double Tmin, double Tmax, const mat & DN, double delta, uvec &low)
{
///
/// \fn mat computeG_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec &low)
/// \brief Returns the matrix G of size (1+M*K)-by-(1+M*K)
/// \param M number of neurons
/// \param  K number of windows
/// \param  Tmin beginning time
/// \param  Tmax end time
/// \param  DN matrix DataNeur, matrix of size M-by-(nmax+1), nmax is the maximum of spikes per neuron
/// \param  delta  K*delta is equal to the scope
/// \param  low array of size Ntot(= number of columns of DS) returned by computeb_V2 function, contains the first index where  
/// \return matrix G of size (1+M*K)-by-(1+M*K)
///
  DataSpike DS(DN);
  mat G = computeG_V2(M, K, Tmin, Tmax, DS, delta, low); 
  return G;
}

