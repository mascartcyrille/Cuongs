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

  int M = 1;
  double Tmin = 0.;
  double Tmax = 2.;
  double delta = 0.02; //0.5;  
  int nmax = 10;
  int K = 5; //10;

  outpath = "./resu3";
  string sim_name = "_3";
  size_t debneur = 1;
  string ficname = "../data/n";
  run_short_test(M, K, delta, Tmin, Tmax, ficname, debneur, outpath, sim_name, nmax);
}
