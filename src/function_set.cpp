#include "function_set.h"
#include "basis_function.h"

namespace dcgp {

/// Append a function to a function set
/**
 *  Append a function to a function set. Duplicate functions are not added.
 */
void function_set::push_back(const basis_function* f)
{
    const typename function_set::size_type size = this->size();

    for(unsigned int i = 0; i < size; i++)
        if(at(i)==f) return;
    std::vector<const basis_function*>::push_back(f);
}

/// Append another function set to a function set
/**
 *  Merge another function set into this one by adding new functions at the end.
 *  Duplicate functions are not added. The order of both function sets is
 *  preserved.
 */
void function_set::push_back(const function_set& f)
{
    const typename function_set::size_type size = f.size();

    reserve(this->size()+size);
    for(unsigned int i = 0; i < size; i++)
        push_back(f.at(i));
}

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

/// Merge two function sets
/**
 *  Merge two function sets by concatenation. Duplicate functions of the second
 *  set already present in the first are not added.
 *  The order of both function sets is preserved.
 */
function_set operator+(const function_set& a, const function_set& b)
{
    function_set res(a);
    res.push_back(b);

    return res;
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
