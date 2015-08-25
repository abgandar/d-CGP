#include <cmath>
#include <string>
#include <DA/DA.h>

#include "wrapped_functions.h"

namespace dcgp {

double my_log_sigmoid(double x, double y)
{
    (void)y;
    return 1.0/(1.0+std::exp(-x));
}

DACE::DA my_d_log_sigmoid(const DACE::DA& x, const DACE::DA& y)
{
    (void)y;
    return 1.0/(1.0+DACE::exp(-x));
}

std::string my_print_log_sigmoid(const std::string& s1, const std::string& s2, bool simplify)
{
    (void)s2;
    (void)simplify;
    return "(1.0/(1.0+exp(-"+s1+")))";
}

} // dcgp namespace ends
