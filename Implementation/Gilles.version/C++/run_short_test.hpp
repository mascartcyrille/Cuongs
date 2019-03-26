#include <chrono>
#include <cstdlib>
#include <sstream>
#include "misc.hpp"
#include "DataSpike.hpp"
#include "computeb_V2.hpp"
#include "computeG_V2.hpp"
#include "computeGk1k2_V2.hpp"

mat run_and_time(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec & low, mat (*f)(int, int, double, double, const DataSpike &, double, uvec &),string varname);

umat run_and_time(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec & low, uvec & cnt, umat (*f)(int, int, double, double, const DataSpike &, double, uvec &, uvec &),string varname);

void run_short_test(int M, int K, double delta, double Tmin, double Tmax, string ficname, size_t debneur, string outpath, string sim_name, size_t nmax=600);
