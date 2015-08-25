#ifndef DCGP_FUNCTION_SET_H
#define DCGP_FUNCTION_SET_H

#include <vector>
#include "basis_function.h"

namespace dcgp {

/// Function set
/**
 * Contains, as static members, several std::vector of dcgp::basis_function containing 
 * function sets of common use. The user can access each set via the syntax 
 * dcgp::function_set::SET_NAME
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class function_set : public std::vector<const basis_function*>
{
public:
	function_set(const function_set& f) : std::vector<const basis_function*>(f) {};
	function_set(std::initializer_list<const basis_function*> l) : std::vector<const basis_function*>(l) {};

    /// Convenience function to return a reference to a basis function directly instead of a pointer
    const basis_function& operator[](const typename function_set::size_type i) const { return *(std::vector<const basis_function*>::operator[](i)); };
    /// Append a function to a function set
    void push_back(const basis_function* f);
    /// Append another function set to a function set
    void push_back(const function_set& f);

    /// set of all predefined basis functions
    static const function_set all;
    /// set of the basic basis functions (+,-,*,/)
    static const function_set basic;
    /// set of the extended basis functions (sqrt,pow,exp,log)
    static const function_set extended;
    /// set of the trigonometric basis functions (sin,cos,tan,asin,acos,atan)
    static const function_set trig;
    /// set of the hyperbolic basis functions (sinh,cosh,tanh)
    static const function_set hyp;
};

std::ostream &operator<<(std::ostream &, const function_set &);
function_set operator+(const function_set& a, const function_set& b);

} // end of namespace dcgp

#endif // DCGP_FUNCTION_SET_H
 