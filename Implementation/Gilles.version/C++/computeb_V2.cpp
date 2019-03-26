#include "misc.hpp"
#include "DataSpike.hpp"

umat computeb_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec & low, uvec & cnt)
{
  const vec & T = DS._T;
  const uvec & neur = DS._neur;
  int Ntot = T.size();
  //  cout << Ntot << endl;
  double A = K*delta;
  //  cout << A << endl;
  umat b((1+M*K),M, fill::zeros);
  low.zeros(Ntot);
  cnt.zeros(M);
  int ilow = 0;
  double eps = 1e-12;
  double t; int r;
  int k, l;
  for(int i=0;i<Ntot;i++)
    {
      t = T(i); r = neur(i);
      //      cout << t << " " << r << endl;
      while(abs(T(ilow)-t)>A)
	ilow++;
      low(i) = ilow;
      //      cout << ilow<< endl;
      if((Tmin<t)&&(t<=Tmax))
	{
	  //	  cout << "verif cnt " <<endl;
	  //	  cout << Tmin << " " << t << " " << Tmax<< endl;
	  cnt(r-1)++;
	  //
	  //	}
	  for(int j=ilow;j<i;j++)
	    {
	      if(abs(t-T(j))<=eps) break;
	      k = get_k(t-T(j), delta);
	      l = neur(j);
	      b((l-1)*K+k,r-1) +=1;
	    }
	}
    }
  for(int i=0;i<M;i++)
    b(0,i) = cnt(i); 
  return b;
}


umat computeb_V2(int M, int K, double Tmin, double Tmax, const mat & DN, double delta, uvec & low, uvec & cnt)
{
  DataSpike DS(DN);
  umat b = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  return b;
}
