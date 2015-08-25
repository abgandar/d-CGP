#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include <DA/DA.h>
#include <boost/lexical_cast.hpp>

#include "exceptions.h"

namespace dcgp {

/// Basis function
/**
 * This struct represent a generic function (or expression) in d-CGP. It contains
 * virtual member functions that allow to compute the function
 * value, its derivatives and its symbolic representation.
 *
 * All functions that are represented in a d-CGP encoding must derive from this class.
 *
 * @author Alexander Wittig (alexander.wittig@esa.int)
 */
class basis_function
{
    public:
    /// type to indicate how many connections a basis_function can make (limited to maximum of two by design)
    enum fn_type { FN_CONST, FN_UNARY, FN_BINARY };

    /// Constructor from name and number of connections
    basis_function(std::string name, fn_type type = FN_BINARY):m_name(name),m_type(type) {}

    /// Overload of operator(double, double)
    /**
    * Allows to call a dcgp::basis_function with the syntax f(double x, double y) and get
    * the function value in x,y in return
    */
    virtual double operator()(double x, double y) const = 0;

    /// Overload of operator(DA,DA)
    /**
    * Allows to call a dcgp::basis_function with the syntax f(DA x, DA y) and get
    * the function value and all derivatives up to the current computation order in x,y in return
    */
    virtual DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const = 0;

    /// Overload of operator(std::string, std::string)
    /**
    * Allows to call a dcgp::basis_function with the syntax f(std::string x, std::string y) and get
    * a symbolic representation of the function in the variables std::string x, std::string y
    */
    virtual std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const = 0;

    /// Its name
    std::string m_name;
    /// Function type (number of input arguments)
    fn_type m_type;
};

std::ostream& operator<<(std::ostream& os, const basis_function& obj);

using my_fun_type = std::function<double(double, double)>;
using d_my_fun_type = std::function<DACE::DA(const DACE::DA&, const DACE::DA&)>;
using my_print_fun_type = std::function<std::string(const std::string, const std::string, bool)>;

/// Basis functor
/**
 * Basic function containing std::function of type my_fun_type, my_d_fun_type,
 * my_print_fun_type, that allow to compute the function value, its derivatives
 * and its symbolic representation.
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Alexander Wittig (alexander.wittig@esa.int)
 */
class basis_functor : public basis_function
{
    public:
    /// Constructor from std::function construction arguments
    template <typename T, typename U, typename V>
    basis_functor(T &&f, U &&df, V&&pf, std::string name, fn_type type = FN_BINARY):basis_function(name,type), m_f(std::forward<T>(f)), m_df(std::forward<U>(df)), m_pf(std::forward<V>(pf)) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;

    /// The function
    my_fun_type m_f;
    /// The DA function
    d_my_fun_type m_df;
    /// Its symbolic representation
    my_print_fun_type m_pf;
};

/// Constant
class basis_cons : public basis_function
{
    public:
    /// Constructor
    basis_cons(const double val = 1.0):basis_function("cons("+boost::lexical_cast<std::string>(val)+")",FN_CONST), m_val(val) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;

    double m_val;
};
/// Global instance of the zero constant
static const basis_cons zero(0.0);
/// Global instance of the one constant
static const basis_cons one(1.0);

/// Addition
class basis_sum : public basis_function
{
    public:
    /// Constructor
    basis_sum():basis_function("sum",FN_BINARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the addition
static const basis_sum sum;

/// Subtraction
class basis_diff : public basis_function
{
    public:
    /// Constructor
    basis_diff():basis_function("diff",FN_BINARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the subtraction
static const basis_diff diff;

/// Multiplication
class basis_mul : public basis_function
{
    public:
    /// Constructor
    basis_mul():basis_function("mul",FN_BINARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the multiplication
static const basis_mul mul;

/// Division
class basis_div : public basis_function
{
    public:
    /// Constructor
    basis_div():basis_function("div",FN_BINARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the division
static const basis_div div;

/// Power
class basis_pow : public basis_function
{
    public:
    /// Constructor
    basis_pow():basis_function("pow",FN_BINARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the power
static const basis_pow pow;

/// Square root
class basis_sqrt : public basis_function
{
    public:
    /// Constructor
    basis_sqrt():basis_function("sqrt",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the square root
static const basis_sqrt sqrt;

/// Exponential
class basis_exp : public basis_function
{
    public:
    /// Constructor
    basis_exp():basis_function("exp",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the exponential
static const basis_exp exp;

/// Natural logarithm
class basis_log : public basis_function
{
    public:
    /// Constructor
    basis_log():basis_function("log",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the natural logarithm
static const basis_log log;


/// Sine
class basis_sin : public basis_function
{
    public:
    /// Constructor
    basis_sin():basis_function("sin",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the sine
static const basis_sin sin;

/// Cosine
class basis_cos : public basis_function
{
    public:
    /// Constructor
    basis_cos():basis_function("cos",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the cosine
static const basis_cos cos;

/// Tangent
class basis_tan : public basis_function
{
    public:
    /// Constructor
    basis_tan():basis_function("tan",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the tangent
static const basis_tan tan;

/// Arcsine
class basis_asin : public basis_function
{
    public:
    /// Constructor
    basis_asin():basis_function("asin",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the arcsine
static const basis_asin asin;

/// Arccosine
class basis_acos : public basis_function
{
    public:
    /// Constructor
    basis_acos():basis_function("acos",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the arccosine
static const basis_acos acos;

/// Arctan
class basis_atan : public basis_function
{
    public:
    /// Constructor
    basis_atan():basis_function("atan",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the arctangent
static const basis_atan atan;

/// Hyperbolic sine
class basis_sinh : public basis_function
{
    public:
    /// Constructor
    basis_sinh():basis_function("sinh",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the hyperbolic sine
static const basis_sinh sinh;

/// Hyperbolic cosine
class basis_cosh : public basis_function
{
    public:
    /// Constructor
    basis_cosh():basis_function("cosh",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the hyperbolic cosine
static const basis_cosh cosh;

/// Hyperbolic tangent
class basis_tanh : public basis_function
{
    public:
    /// Constructor
    basis_tanh():basis_function("tanh",FN_UNARY) {}

    /// Overload of operator(double, double)
    double operator()(double x, double y) const;

    /// Overload of operator(DA,DA)
    DACE::DA operator()(const DACE::DA& x, const DACE::DA& y) const;

    /// Overload of operator(std::string, std::string)
    std::string operator()(const std::string s1, const std::string s2, bool simplify = false) const;
};
/// Global instance of the hyperbolic tangent
static const basis_tanh tanh;

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
 