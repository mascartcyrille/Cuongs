#ifndef __COMPUTEGK1K2_V2_HPP__
#define __COMPUTEGK1K2_V2_HPP__
#include "computeb_V2.hpp"
mat computeGk1k2_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec &low);
mat computeGk1k2_V2(int M, int K, double Tmin, double Tmax, const mat & DN, double delta, uvec &low);
#endif
