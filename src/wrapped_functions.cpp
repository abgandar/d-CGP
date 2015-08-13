#include <cmath>
#include <string>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>

#include "wrapped_functions.h"
#include "std_overloads.h"


namespace dcgp {

double my_sum(double x, double y)
{
        return x + y;
}

DACE::DA d_my_sum(const DACE::DA& x, const DACE::DA& y)
{
    return x + y;
}

std::string print_my_sum(const std::string& s1, const std::string& s2)
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
    return ("(" + s1 + "+" + s2 + ")");
}


double my_diff(double x, double y)
{
        return x - y;
}

DACE::DA d_my_diff(const DACE::DA& x, const DACE::DA& y)
{
    return x - y;
}

std::string print_my_diff(const std::string& s1, const std::string& s2)
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
    return ("(" + s1 + "-" + s2 + ")");
}

double my_mul(double x, double y)
{
        return (x * y);
}

DACE::DA d_my_mul(const DACE::DA& x, const DACE::DA& y)
{
    return x * y;
}

std::string print_my_mul(const std::string& s1, const std::string& s2)
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
    return ("(" + s1 + "*" + s2 + ")");
}

double my_div(double x, double y)
{
        return x / y;
}

DACE::DA d_my_div(const DACE::DA& x, const DACE::DA& y)
{
    if (cons(y)==0.0)
    {
        throw derivative_error("Derivative of div does not exist at this point");
    }
    return x / y;
}

std::string print_my_div(const std::string& s1, const std::string& s2)
{
    if (s1 == "0" && s2 != "0")
    {
        return "0";
    }
    else if (s1 == s2 && s1 != "0")
    {
        return "1";
    }
    return ("(" + s1 + "/" + s2 + ")");
}

double my_pow(double b, double c)
{
        return pow(fabs(b),c);
}

DACE::DA d_my_pow(const DACE::DA& x, const DACE::DA& y)
{
    // We derive this by setting a = exp(c * ln(|b|))
    // TODO (AW): maybe an arbitrary power is not the best function to include. Instead maybe including exp() and ln() is sufficient?
    if (cons(x)<0)
    {
        return exp(y*log(-x));
    }
    else if (cons(x)>0)
    {
        return exp(y*log(x));
    }
    else
    {
        throw derivative_error("Derivative of pow does not exist at this point");
    }
}

std::string print_my_pow(const std::string& s1, const std::string& s2)
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
    return ("abs(" + s1 + ")^(" + s2 + ")");
}

double my_sqrt(double b, double c)
{
        (void)c;
        return sqrt(fabs(b));
}

DACE::DA d_my_sqrt(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    if (cons(b)<0)
    {
        return sqrt(-b);
    }
    else if (cons(b)>0)
    {
        return sqrt(b);
    }
    else
    {
        throw derivative_error("Derivative of sqrt does not exist at this point");
    }
}

std::string print_my_sqrt(const std::string& s1, const std::string& s2)
{
    (void)s2;
    if (s1 == "0")
    {
        return "0";
    }
    else if (s1 == "1")
    {
        return "1";
    }
    return ("sqrt(abs(" + s1 + "))");
}

double my_exp(double b, double c)
{
        (void)c;
        return exp(b);
}

DACE::DA d_my_exp(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return exp(b);
}

std::string print_my_exp(const std::string& s1, const std::string& s2)
{
    (void)s2;
    if (s1 == "0")
    {
        return "1";
    }
    return ("exp(" + s1 + ")");
}

double my_log(double b, double c)
{
        (void)c;
        return log(fabs(b));
}

DACE::DA d_my_log(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    if (cons(b)<0)
    {
        return log(-b);
    }
    else if (cons(b)>0)
    {
        return log(b);
    }
    else
    {
        throw derivative_error("Derivative of log does not exist at this point");
    }
}

std::string print_my_log(const std::string& s1, const std::string& s2)
{
    (void)s2;
    if (s1 == "1")
    {
        return "0";
    }
    return ("log(abs(" + s1 + "))");
}

double my_sin(double b, double c)
{
        (void)c;
        return sin(b);
}

DACE::DA d_my_sin(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return sin(b);
}

std::string print_my_sin(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("sin(" + s1 + ")");
}

double my_cos(double b, double c)
{
        (void)c;
        return cos(b);
}

DACE::DA d_my_cos(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return cos(b);
}

std::string print_my_cos(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("cos(" + s1 + ")");
}

double my_tan(double b, double c)
{
        (void)c;
        return tan(b);
}

DACE::DA d_my_tan(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return tan(b);
}

std::string print_my_tan(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("tan(" + s1 + ")");
}

double my_asin(double b, double c)
{
        (void)c;
        return asin(b);
}

DACE::DA d_my_asin(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    if (fabs(cons(b))>=1.0)
    {
        throw derivative_error("Derivative of asin does not exist at this point");
    }
    return asin(b);
}

std::string print_my_asin(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("asin(" + s1 + ")");
}

double my_acos(double b, double c)
{
        (void)c;
        return acos(b);
}

DACE::DA d_my_acos(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    if (fabs(cons(b))>=1.0)
    {
        throw derivative_error("Derivative of acos does not exist at this point");
    }
    return acos(b);
}

std::string print_my_acos(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("acos(" + s1 + ")");
}

double my_atan(double b, double c)
{
        (void)c;
        return atan(b);
}

DACE::DA d_my_atan(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return atan(b);
}

std::string print_my_atan(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("atan(" + s1 + ")");
}

double my_sinh(double b, double c)
{
        (void)c;
        return sinh(b);
}

DACE::DA d_my_sinh(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return sinh(b);
}

std::string print_my_sinh(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("sinh(" + s1 + ")");
}

double my_cosh(double b, double c)
{
        (void)c;
        return cosh(b);
}

DACE::DA d_my_cosh(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return cosh(b);
}

std::string print_my_cosh(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("cosh(" + s1 + ")");
}

double my_tanh(double b, double c)
{
        (void)c;
        return tanh(b);
}

DACE::DA d_my_tanh(const DACE::DA& b, const DACE::DA& c)
{
    (void)c;
    return tanh(b);
}

std::string print_my_tanh(const std::string& s1, const std::string& s2)
{
    (void)s2;
    return ("tanh(" + s1 + ")");
}

DACE::DA d_not_implemented(const DACE::DA& b, const DACE::DA& c)
{
    (void)b;
    (void)c;
    throw derivative_error("Differentiation has not been implemented ... you can use CGP but not d-CGP");
}

} // dcgp namespace ends
