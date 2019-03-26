#ifndef __MISC_HPP__
#define __MISC_HPP__

#include <armadillo>
#include <cmath>
#include <iostream>
using namespace::std;
using namespace::arma;

int get_k(double x, double delta, double eps=1e-12);
unsigned int get_low_index(double val, const vec& T);
unsigned int get_up_index(double val, const vec& T);
umat get_k1k2(size_t K, double x, double delta, double eps=1e-12);
size_t get_digit(size_t n);
#endif
