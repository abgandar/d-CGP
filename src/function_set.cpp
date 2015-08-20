#include "function_set.h"
#include "basis_function.h"

namespace dcgp {

/// Overload stream insertion operator for function_set
std::ostream &operator<<(std::ostream &os, const function_set &fs)
{
    const typename function_set::size_type len = fs.size();
    os << '[';
    for (typename function_set::size_type i = 0; i < len; ++i) {
        os << fs[i].m_name;
        if (i != len-1) {
            os << ", ";
        }
    }
    os << ']';

    return os;
}

/// set of all predefined basis functions
const function_set function_set::all({&zero, &one, &sum, &diff, &mul, &div, &sqrt, &pow, &exp, &log, &sin, &cos, &tan, &asin, &acos, &atan, &sinh, &cosh, &tanh});
/// set of the basic basis functions (+,-,*,/)
const function_set function_set::basic({&zero, &one, &sum, &diff, &mul, &div});
/// set of the extended basis functions (sqrt,pow,exp,log)
const function_set function_set::extended({&sqrt, &pow, &exp, &log});
/// set of the trigonometric basis functions (sin,cos,tan,asin,acos,atan)
const function_set function_set::trig({&sin, &cos, &tan, &asin, &acos, &atan});
/// set of the hyperbolic basis functions (sinh,cosh,tanh)
const function_set function_set::hyp({&sinh, &cosh, &tanh});

} // end of namespace dcgp
