#include "function_set.h"
#include "wrapped_functions.h"
#include "exceptions.h"

namespace dcgp {

function_set::function_set() : m_functions() {};

function_set::function_set(const std::vector<std::string>& list)
{
    for (auto function_name : list)
    {
        push_back(function_name);
    }
}

void function_set::push_back(const std::string& function_name)
{
    if (function_name=="sum")
        m_functions.emplace_back(my_sum,d_my_sum,print_my_sum, function_name);
    else if (function_name=="diff")
        m_functions.emplace_back(my_diff,d_my_diff,print_my_diff, function_name);
    else if (function_name=="mul")
        m_functions.emplace_back(my_mul,d_my_mul,print_my_mul, function_name);
    else if (function_name=="div")
        m_functions.emplace_back(my_div,d_my_div,print_my_div, function_name);
    else if (function_name=="sqrt")
        m_functions.emplace_back(my_sqrt,d_my_sqrt,print_my_sqrt, function_name);
    else if (function_name=="pow")
        m_functions.emplace_back(my_pow,d_my_pow,print_my_pow, function_name);
    else if (function_name=="exp")
        m_functions.emplace_back(my_exp,d_my_exp,print_my_exp, function_name);
    else if (function_name=="log")
        m_functions.emplace_back(my_log,d_my_log,print_my_log, function_name);
    else if (function_name=="sin")
        m_functions.emplace_back(my_sin,d_my_sin,print_my_sin, function_name);
    else if (function_name=="cos")
        m_functions.emplace_back(my_cos,d_my_cos,print_my_cos, function_name);
    else if (function_name=="tan")
        m_functions.emplace_back(my_tan,d_my_tan,print_my_tan, function_name);
    else if (function_name=="asin")
        m_functions.emplace_back(my_asin,d_my_asin,print_my_asin, function_name);
    else if (function_name=="acos")
        m_functions.emplace_back(my_acos,d_my_acos,print_my_acos, function_name);
    else if (function_name=="atan")
        m_functions.emplace_back(my_atan,d_my_atan,print_my_atan, function_name);
    else if (function_name=="sinh")
        m_functions.emplace_back(my_sinh,d_my_sinh,print_my_sinh, function_name);
    else if (function_name=="cosh")
        m_functions.emplace_back(my_cosh,d_my_cosh,print_my_cosh, function_name);
    else if (function_name=="tanh")
        m_functions.emplace_back(my_tanh,d_my_tanh,print_my_tanh, function_name);
    else 
        throw input_error("Unimplemented function " + function_name);
}

void function_set::clear()
{
    m_functions.clear();
}

std::vector<basis_function> function_set::operator()() const
{
    return m_functions;
}


} // end of namespace dcgp
