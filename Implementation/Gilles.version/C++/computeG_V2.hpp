#ifndef __COMPUTEG_V2_HPP__
#define __COMPUTEG_V2_HPP__
mat computeG_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec &low);
mat computeG_V2(int M, int K, double Tmin, double Tmax, const mat & DN, double delta, uvec &low);
#endif
