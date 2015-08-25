#include <string>
#include <iostream>

#include "basis_function.h"

namespace dcgp {

/// Overload stream operator for dcgp::basis_function
/**
 * @param[in] obj dcgp::basis_function to be inserted into the stream.
 * @param[out] os std::ostream to which the problem will be streamed.
 *
 * @return reference to os.
 */
std::ostream& operator<<(std::ostream& os, const basis_function& obj)
{
    os << obj.m_name;
    return os;
}

/// Overload of operator(double, double)
/**
* Call the associated functor for function evaluation at x,y
*/
double basis_functor::operator()(double x, double y) const
{
    return m_f(x,y);
}

/// Overload of operator(DA,DA)
/**
* Call the associated functor for DA function (derivative) evaluation at x,y
*/
DACE::DA basis_functor::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    return m_df(x,y);
}

/// Overload of operator(std::string, std::string, bool)
/**
* Call the associated functor for function printing
*/
std::string basis_functor::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    return m_pf(s1,s2,simplify);
}

/// Overload of operator(double, double)
double basis_cons::operator()(double x, double y) const 
{
    (void)x;
    (void)y;
    return m_val;
}

/// Overload of operator(DA,DA)
DACE::DA basis_cons::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void)x;
    (void)y;
    return m_val;
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_cons::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s1;
    (void)s2;
    (void)simplify;
    return boost::lexical_cast<std::string>(m_val);
}

/// Overload of operator(double, double)
double basis_sum::operator()(double x, double y) const 
{
    return x+y;
}

/// Overload of operator(DA,DA)
DACE::DA basis_sum::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    return x+y;
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_sum::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    if(simplify)
    {
        if (s1 == s2) 
        {
            return "(2*"+s1+")";
        }
        else if (s1 == "0")
        {
            return s2;
        }
        else if (s2 == "0")
        {
            return s1;
        }
    }
    return ("(" + s1 + "+" + s2 + ")");
}

/// Overload of operator(double, double)
double basis_diff::operator()(double x, double y) const 
{
    return x-y;
}

/// Overload of operator(DA,DA)
DACE::DA basis_diff::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    return x-y;
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_diff::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    if(simplify)
    {
        if (s1 == s2) 
        {
            return "0";
        }
        else if (s1 == "0")
        {
            return "(-" + s2 + ")";
        }
        else if (s2 == "0")
        {
            return s1;
        }
    }
    return ("(" + s1 + "-" + s2 + ")");
}

/// Overload of operator(double, double)
double basis_mul::operator()(double x, double y) const 
{
    return x*y;
}

/// Overload of operator(DA,DA)
DACE::DA basis_mul::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    return x*y;
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_mul::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    if(simplify)
    {
        if (s1 == "0" || s2 == "0")
        {
            return "0";
        }
        else if (s1 == s2)
        {
            return s1 + "^2";
        }
        else if (s1 == "1")
        {
            return s2;
        }
        else if (s2 == "1")
        {
            return s1;
        }
    }
    return ("(" + s1 + "*" + s2 + ")");
}

/// Overload of operator(double, double)
double basis_div::operator()(double x, double y) const 
{
    return x/y;
}

/// Overload of operator(DA,DA)
DACE::DA basis_div::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    if (cons(y)==0.0)
    {
        throw derivative_error("Derivative of div does not exist at this point");
    }
    return x/y;
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_div::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    if(simplify)
    {
        if (s1 == "0" && s2 != "0")
        {
            return "0";
        }
        else if (s1 == s2 && s1 != "0")
        {
            return "1";
        }
    }
    return ("(" + s1 + "/" + s2 + ")");
}

/// Overload of operator(double, double)
double basis_pow::operator()(double x, double y) const 
{
    return std::pow(std::abs(x),y);
}

/// Overload of operator(DA,DA)
DACE::DA basis_pow::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    if (cons(x)<0)
    {
        return DACE::exp(y*DACE::log(-x));
    }
    else if (cons(x)>0)
    {
        return DACE::exp(y*DACE::log(x));
    }
    else
    {
        throw derivative_error("Derivative of pow does not exist at this point");
    }
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_pow::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    if(simplify)
    {
        if (s1 == "0" && s2 != "0")
        {
            return "0";
        }
        else if (s1 == "1")
        {
            return "1";
        }
        else if (s2 == "0" && s1 != "0")
        {
            return "1";
        }
        else if (s2 == "1")
        {
            return ("abs(" + s1 + ")");
        }
    }
    return ("abs(" + s1 + ")^(" + s2 + ")");
}

/// Overload of operator(double, double)
double basis_sqrt::operator()(double x, double y) const 
{
    (void)y;
    return std::sqrt(std::abs(x));
}

/// Overload of operator(DA,DA)
DACE::DA basis_sqrt::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void)y;
    if (cons(x)<0)
    {
        return DACE::sqrt(-x);
    }
    else if (cons(x)>0)
    {
        return DACE::sqrt(x);
    }
    else
    {
        throw derivative_error("Derivative of sqrt does not exist at this point");
    }
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_sqrt::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    if(simplify)
    {
        if (s1 == "0")
        {
            return "0";
        }
        else if (s1 == "1")
        {
            return "1";
        }
    }
    return ("sqrt(abs(" + s1 + "))");
}

/// Overload of operator(double, double)
double basis_exp::operator()(double x, double y) const 
{
    (void) y;
    return std::exp(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_exp::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::exp(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_exp::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    if(simplify)
    {
        if (s1 == "0")
        {
            return "1";
        }
    }
    return ("exp(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_log::operator()(double x, double y) const 
{
    (void) y;
    return std::log(std::abs(x));
}

/// Overload of operator(DA,DA)
DACE::DA basis_log::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    if (cons(x)<0)
    {
        return DACE::log(-x);
    }
    else if (cons(x)>0)
    {
        return DACE::log(x);
    }
    else
    {
        throw derivative_error("Derivative of log does not exist at this point");
    }
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_log::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    if(simplify)
    {
        if (s1 == "1")
        {
            return "0";
        }
    }
    return ("log(abs(" + s1 + "))");
}

/// Overload of operator(double, double)
double basis_sin::operator()(double x, double y) const 
{
    (void) y;
    return std::sin(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_sin::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::sin(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_sin::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("sin(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_cos::operator()(double x, double y) const 
{
    (void) y;
    return std::cos(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_cos::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::cos(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_cos::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("cos(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_tan::operator()(double x, double y) const 
{
    (void) y;
    return std::tan(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_tan::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::tan(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_tan::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("tan(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_asin::operator()(double x, double y) const 
{
    (void) y;
    return std::asin(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_asin::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void)y;
    if (std::abs(cons(x))>=1.0)
    {
        throw derivative_error("Derivative of asin does not exist at this point");
    }
    return DACE::asin(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_asin::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("asin(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_acos::operator()(double x, double y) const 
{
    (void) y;
    return std::acos(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_acos::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void)y;
    if (std::abs(cons(x))>=1.0)
    {
        throw derivative_error("Derivative of acos does not exist at this point");
    }
    return DACE::acos(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_acos::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("acos(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_atan::operator()(double x, double y) const 
{
    (void) y;
    return std::atan(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_atan::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::atan(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_atan::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("atan(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_sinh::operator()(double x, double y) const 
{
    (void) y;
    return std::sinh(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_sinh::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::sinh(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_sinh::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("sinh(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_cosh::operator()(double x, double y) const 
{
    (void) y;
    return std::cosh(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_cosh::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::cosh(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_cosh::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("cosh(" + s1 + ")");
}

/// Overload of operator(double, double)
double basis_tanh::operator()(double x, double y) const 
{
    (void) y;
    return std::tanh(x);
}

/// Overload of operator(DA,DA)
DACE::DA basis_tanh::operator()(const DACE::DA& x, const DACE::DA& y) const
{
    (void) y;
    return DACE::tanh(x);
}

/// Overload of operator(std::string, std::string, bool)
std::string basis_tanh::operator()(const std::string s1, const std::string s2, bool simplify) const
{
    (void)s2;
    (void)simplify;
    return ("tanh(" + s1 + ")");
}

} // end of namespace dcgp

 