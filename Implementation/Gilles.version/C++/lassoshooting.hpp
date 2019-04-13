#ifndef __LASSOSHOOTING_HPP__
#define __LASSOSHOOTING_HPP__

#include <armadillo>
#include <iostream>
#include <iomanip>
using namespace::std;
using namespace::arma;

mat lasso_shooting(const mat &G, const umat & b, const mat & d, const vec & a0, int nitmax, double eps, uvec & nit);

//void lasso_shooting1(const mat &G, const uvec & b, const mat & d, const vec & a0, int nitmax, double eps, vec & a, size_t &  nit);
void lasso_shooting1(const mat &G, const uvec & b, const mat & d, const vec & a0, int nitmax, double eps, vec & a, long long unsigned int *  nit);
//void lasso_shooting1(const mat &G, const uvec & b, const mat & d, const vec & a0, int nitmax, double eps, vec & a, size_t nit);

double J(const mat &G, const uvec & b, const uvec & d, const vec & a);
#endif
