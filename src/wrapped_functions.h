#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <string>
#include <vector>

#include <DACE/DA.h>

#include "exceptions.h"

namespace dcgp {

/*--------------------------------------------------------------------------
*                                  BINARY FUNCTIONS
*------------------------------------------------------------------------**/
// f = b + c
double my_sum(double b, double c);
DACE::DA d_my_sum(const DACE::DA& b, const DACE::DA& c);
std::string print_my_sum(const std::string& s1, const std::string& s2);

// f = b - c
double my_diff(double b, double c);
DACE::DA d_my_diff(const DACE::DA& b, const DACE::DA& c);
std::string print_my_diff(const std::string& s1, const std::string& s2);

// f = b * c
double my_mul(double b, double c);
DACE::DA d_my_mul(const DACE::DA& b, const DACE::DA& c);
std::string print_my_mul(const std::string& s1, const std::string& s2);

// f = b / c
double my_div(double b, double c);
DACE::DA d_my_div(const DACE::DA& b, const DACE::DA& c);
std::string print_my_div(const std::string& s1, const std::string& s2);

// f = pow(|b|,c)
double my_pow(double b, double c);
DACE::DA d_my_pow(const DACE::DA& b, const DACE::DA& c);
std::string print_my_pow(const std::string& s1, const std::string& s2);

/*--------------------------------------------------------------------------
*                                  UNARY FUNCTIONS
*------------------------------------------------------------------------**/

// f = sqrt(|b|)
double my_sqrt(double b, double c);
DACE::DA d_my_sqrt(const DACE::DA& b, const DACE::DA& c);
std::string print_my_sqrt(const std::string& s1, const std::string& s2);

// f = exp(b)
double my_exp(double b, double c);
DACE::DA d_my_exp(const DACE::DA& b, const DACE::DA& c);
std::string print_my_exp(const std::string& s1, const std::string& s2);

// f = log(|b|)
double my_log(double b, double c);
DACE::DA d_my_log(const DACE::DA& b, const DACE::DA& c);
std::string print_my_log(const std::string& s1, const std::string& s2);

// f = sin(b)
double my_sin(double b, double c);
DACE::DA d_my_sin(const DACE::DA& b, const DACE::DA& c);
std::string print_my_sin(const std::string& s1, const std::string& s2);

// f = cos(b)
double my_cos(double b, double c);
DACE::DA d_my_cos(const DACE::DA& b, const DACE::DA& c);
std::string print_my_cos(const std::string& s1, const std::string& s2);

// f = tan(b)
double my_tan(double b, double c);
DACE::DA d_my_tan(const DACE::DA& b, const DACE::DA& c);
std::string print_my_tan(const std::string& s1, const std::string& s2);

// f = asin(b)
double my_asin(double b, double c);
DACE::DA d_my_asin(const DACE::DA& b, const DACE::DA& c);
std::string print_my_asin(const std::string& s1, const std::string& s2);

// f = acos(b)
double my_acos(double b, double c);
DACE::DA d_my_acos(const DACE::DA& b, const DACE::DA& c);
std::string print_my_acos(const std::string& s1, const std::string& s2);

// f = atan(b)
double my_atan(double b, double c);
DACE::DA d_my_atan(const DACE::DA& b, const DACE::DA& c);
std::string print_my_atan(const std::string& s1, const std::string& s2);

// f = sinh(b)
double my_sinh(double b, double c);
DACE::DA d_my_sinh(const DACE::DA& b, const DACE::DA& c);
std::string print_my_sinh(const std::string& s1, const std::string& s2);

// f = cosh(b)
double my_cosh(double b, double c);
DACE::DA d_my_cosh(const DACE::DA& b, const DACE::DA& c);
std::string print_my_cosh(const std::string& s1, const std::string& s2);

// f = tanh(b)
double my_tanh(double b, double c);
DACE::DA d_my_tanh(const DACE::DA& b, const DACE::DA& c);
std::string print_my_tanh(const std::string& s1, const std::string& s2);

/*--------------------------------------------------------------------------
*                                  HELPER FUNCTIONS
*------------------------------------------------------------------------**/

DACE::DA d_not_implemented(const DACE::DA& b, const DACE::DA& c);

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H