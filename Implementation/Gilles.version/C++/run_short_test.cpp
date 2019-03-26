#include <chrono>
#include <cstdlib>
#include "misc.hpp"
#include "DataSpike.hpp"
#include "computeb_V2.hpp"
#include "computeG_V2.hpp"
#include "computeGk1k2_V2.hpp"

// for computeG functions (computeG_V2, computeGpw_V2, computeGk1k2_V2)
//
mat run_and_time(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec & low, mat (*f)(int, int, double, double, const DataSpike &, double, uvec &),string varname)
{
  auto start = chrono::high_resolution_clock::now();
  mat G = (*f)(M, K, Tmin, Tmax, DS, delta, low);
  auto finish = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = finish - start;

  cout << "Elapsed time to compute " << varname << ":" << elapsed.count() << " s\n";
  return G;
}

// for computeb_V2 function
//
umat run_and_time(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec & low, uvec & cnt, umat (*f)(int, int, double, double, const DataSpike &, double, uvec &, uvec &),string varname)
{
  auto start = chrono::high_resolution_clock::now();
  umat b = (*f)(M, K, Tmin, Tmax, DS, delta, low, cnt);
  auto finish = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = finish - start;

  cout << "Elapsed time to compute " << varname << ":" << elapsed.count() << " s\n";
  return b;
}


void run_short_test(int M, int K, double delta, double Tmin, double Tmax, string ficname, size_t debneur, string outpath, string sim_name, size_t nmax)
{

  cout << endl;
  cout << "Test " << sim_name << " - " << M << " neurons - computation of b and G with V2 method with C++ armadillo "<< endl;
  cout << endl;
  cout << endl;
  uvec low, cnt;
  string varname;
  //  cout << "ficname " << ficname << " M" << M << " debneur " << debneur << " nmax " << nmax << endl;

  outpath.append("/");
  string s1 = "mkdir -p " + outpath;
  system(s1.c_str());
  
  system("date; hostname; id -a; g++ --version|head -1");
  cout << endl;
  cout << "K= " << K << endl;
  DataSpike DS(ficname, M, debneur, nmax);
  
  cout << "Total number of spikes independently of neurons " << DS._T.size() << endl;

  // b computed with computeb_V2 function
  varname = "b";
  umat b1 = run_and_time(M, K, Tmin, Tmax, DS, delta, low, cnt, &computeb_V2, varname);
  ficname = outpath + varname + sim_name + ".txt";
  b1.save(ficname, raw_ascii);
  
  // G computed with computeG_V2 function
  varname = "G";
  mat G1 = run_and_time(M, K, Tmin, Tmax, DS, delta, low, &computeG_V2, varname);
  ficname = outpath + varname + sim_name + ".txt";
  G1.save(ficname, raw_ascii);
  
  // Gk computed with computeGk1k2_V2 function
  varname = "Gk";
  mat Gk = run_and_time(M, K, Tmin, Tmax, DS, delta, low, &computeGk1k2_V2, varname);
  ficname = outpath + varname + sim_name + ".txt";
  Gk.save(ficname, raw_ascii);
}



