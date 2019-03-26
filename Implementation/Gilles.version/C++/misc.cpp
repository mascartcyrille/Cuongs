#include "misc.hpp"
int get_k(double x, double delta, double eps)//=1e-12)
{
//     Computation of k such that x belongs to ](k-1)*delta, k*delta]
  int k = x/delta + 1;
  if(abs(k*delta -x)<eps) 
    return k;
  if(abs((k-1)*delta -x)<eps)
    return k-1;
  return k;
 }

unsigned int get_low_index(double x, const vec& T)
{
  // Computation of ind, the first index in T such that T[ind]>x. T is supposed to be sorted in ascending order.
  // Possible values of ind are between 0 and T.size() (included)
  // if get_low_index==0, it means  that all values of T satisfy T>x
  // if get_low_index==T.size(), it means  that all values of T satisfy T<=x
  int ind = 0;
  int n = T.size();
  while(ind < n) {
    if(T(ind)<=x)
      ind++;
    else
      break;
  }
  return ind;
}

unsigned int get_up_index(double x, const vec& T)
{
  // Computation of ind, the last index in T such that T[ind]<=x """
  // Possible values of ind are between -1 and (T.size()-1) (included)
  int ind = T.size() - 1;
  while(ind >= 0) {
    if(T(ind)>x)
      ind--;
    else
      break;
  }
  return ind;
}

umat get_k1k2(size_t K, double x, double delta, double eps)
{
  // 1st column of Mk contains k1, 2nd col contains k2
  // k2 - k1 = chi
  // k2 - k1 = chi + 1
  // k1 - k2 = chi + 1
  // k1 - k2 = chi
  // we have 1 <= k1, k2 <= K
  // if chi==0, k2 - k1 =  chi gives the same indices as k1 - k2 = chi
  //
  size_t chi = int(x/delta);
  if(fabs((chi+1)*delta - x)< eps)
    chi++;
  //  cout << "get_k1k2 "<< x/delta << " " << chi<<endl;
  umat Mk(K*K, 2);
  size_t irow=0;
  for(size_t n=chi+1;n<K+1;n++)
    {
      // k1-k2 = chi
      Mk(irow,0) = n;
      Mk(irow,1) = n - chi;
      irow++;
      // k2-k1 = chi
      if(chi>0) {
      Mk(irow,0) = n-chi; 
      Mk(irow,1) = n; 
      irow++;
      }
    }
    for(size_t n=chi+2;n<K+1;n++)
    {
      // k1-k2 = chi+1
      Mk(irow,0) = n;
      Mk(irow,1) = n - chi - 1;
      irow++;
      // k2-k1 = chi+1
      Mk(irow,0) = n-chi-1; 
      Mk(irow,1) = n; 
      irow++;
    }
    //    cout << "irow= " << irow << endl;
    if(irow>0)
      return Mk.submat(0, 0, irow-1, 1);
    else
      {
	//	cout << "get_k1k2 "<< x/delta << " " << chi<<" " << K <<endl;
	return zeros<umat>(0, 0);
      }
      }


// returns the numbers of digits of a given positive number
// example: n=2 => returns 1
//          n=25 => returns 2
//          n=525 => returns 3
//          n=0 => returns 1
size_t get_digit(size_t n)
{
  size_t num = 1;
  while((n/10)>0) {
    n = n/10;
    num++;
  }
  return num;
}
