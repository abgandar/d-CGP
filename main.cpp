#include <iostream>
#include <iomanip>
#include <ctime>

#include <DA/DA.h>

#include "src/dcgp.h"

int main() {
    // We define the set of functions we want to use
    dcgp::function_set basic_set = dcgp::function_set::basic;

    // We instantiate a d-CGP expression
    unsigned int n_inputs = 3;
    unsigned int n_outputs = 2;
    unsigned int n_rows = 3;
    unsigned int n_columns = 7;
    unsigned int n_level_backs = 4;
    dcgp::expression simple(n_inputs,n_outputs,n_rows,n_columns,n_level_backs,basic_set);

    // initialize DACE to handle derivatives up to order n_derivatives in n_inputs variables
    unsigned int n_derivatives = 3;
    DACE::DA::init(n_derivatives,n_inputs);

    // We inspect it
    std::cout << simple << std::endl;

    // We compute the expression value and its derivatives in a point
    std::vector<double> in_num({2.,3.,4.});
    std::cout << "Point is:" << in_num << std::endl;
    std::cout << "Numerical value = " << simple(in_num) << std::endl;

    std::vector<DACE::DA> deriv = simple.differentiate(in_num);
    
    std::vector<double> jet_0 = simple.differentiate({0,0},deriv);
    std::cout << "Numerical values d^2/dx^2 = " << jet_0 << std::endl;

    std::vector<double> jet_1 = simple.differentiate({1,1},deriv);
    std::cout << "Numerical values d^2/dy^2 = " << jet_1 << std::endl;

    std::vector<double> jet_2 = simple.differentiate({2,2},deriv);
    std::cout << "Numerical values d^2/dz^2 = " << jet_2 << std::endl;

    std::vector<double> jet_3 = simple.differentiate({0,1},deriv);
    std::cout << "Numerical values d^2/dxdy = " << jet_3 << std::endl;

    std::vector<double> jet_4 = simple.differentiate({1,2},deriv);
    std::cout << "Numerical values d^2/dydz = " << jet_4 << std::endl;

    std::vector<double> jet_5 = simple.differentiate({0,2},deriv);
    std::cout << "Numerical values d^2/dxdz = " << jet_5 << std::endl;

    // We stream a symbolic representation of the expression
    std::vector<std::string> in_sym({"x","y","z"});
    std::cout << "Symbolic value = " << simple(in_sym) << std::endl;

    return 0;
}

/* Possible output:

d-CGP Expression:
    Number of inputs:       3
    Number of outputs:      2
    Number of rows:         2
    Number of columns:      7
    Number of levels-back allowed:  4

    Resulting lower bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  ... ]
    Resulting upper bounds: [3, 2, 2, 3, 2, 2, 3, 4, 4, 3, 4, 4, 3, 6, 6, 3, 6, 6, 3, 8,  ... ]

    Current expression (encoded):   [1, 2, 0, 3, 1, 1, 0, 1, 3, 2, 0, 3, 2, 6, 2, 0, 3, 4, 1, 6,  ... ]
    Active nodes:           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13]
    Active genes:           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  ... ]

Point is:[2, 3, 4]
Numerical value = [0.20000000000000001, 48]
Numerical values d^n/dx^n = [[0.20000000000000001, 48], [0.23999999999999999, -16], [-0.30399999999999999, -24]]
Numerical values d^n/dy^n = [[0.20000000000000001, 48], [-0.040000000000000001, 0], [0.016, 0]]
Numerical values d^n/dz^n = [[0.20000000000000001, 48], [0.16, 52], [-0.064000000000000001, 36]]
Symbolic value = [(((x*(z-x))-((z-x)+1))/(y+(z-x))), (((x*(z-x))*z)*((z-x)+1))]

**/
