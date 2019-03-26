#include "run_short_test.hpp"

int main(int argc, char *argv[])
{
  string outpath;
  if(argc==1)
    {
    outpath = "/workdir/gscarella/lassohawkes/C++/";
    cout << "b and G matrices are written in "<< outpath <<" folder. If it is not ok for you, please add the correct path as an argument of the main function"<<endl;
    }
    else
    outpath = argv[1];

  int M = 100;
  double Tmin = 0.;
  double Tmax = 8.;
  double delta = 0.02; //0.5;  
  int nmax = 600;
  int K = 5; //10;

  string sim_name = "_11";
  size_t debneur = 500;
  string ficname = "../data/M=1000/N";
  run_short_test(M, K, delta, Tmin, Tmax, ficname, debneur, outpath, sim_name, nmax);
}
