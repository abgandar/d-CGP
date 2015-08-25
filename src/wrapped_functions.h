#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <string>
#include <DA/DA.h>

#include "basis_function.h"

namespace dcgp {

/*--------------------------------------------------------------------------
*                        Example wrapped function
*------------------------------------------------------------------------**/
// the logistic sigmoid function (just for demonstration purposes)
double my_log_sigmoid(double b, double c);
DACE::DA my_d_log_sigmoid(const DACE::DA& b, const DACE::DA& c);
std::string my_print_log_sigmoid(const std::string& s1, const std::string& s2, bool simplify = false);

// create global instance of the logistic sigmoid
static const basis_functor log_sigmoid(my_log_sigmoid, my_d_log_sigmoid, my_print_log_sigmoid, "logistic sigmoid", basis_function::FN_UNARY);

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H