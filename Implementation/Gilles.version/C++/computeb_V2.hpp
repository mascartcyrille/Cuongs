#ifndef __COMPUTEB_V2_HPP__
#define __COMPUTEB_V2_HPP__
umat computeb_V2(int M, int K, double Tmin, double Tmax, const DataSpike & DS, double delta, uvec & low, uvec & cnt);
umat computeb_V2(int M, int K, double Tmin, double Tmax, const mat & DN, double delta, uvec & low, uvec & cnt);
#endif
