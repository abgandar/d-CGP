#include <iostream>
#include <sstream>
#include <random>
#include <limits>
#include <cmath>

#include "expression.h"
#include "std_overloads.h"


namespace dcgp {

/// Constructor
/** Constructs a d-cgp expression
 *
 * \param[in] n number of inputs (independent variables)
 * \param[in] m number of outputs (dependent variables)
 * \param[in] r number of rows of the cartesian cgp
 * \param[in] c number of columns of the cartesian cgp
 * \param[in] l number of levels-back allowed for the cartesian cgp
 * \param[in] f function set. An std::vector of dcgp::basis_function
 * \param[in] tol tolerance to be used in case dcgp::expression::HITS_BASED is used as fitness evaluation 
 * \param[in] seed seed for the random number generator (initial expression  and mutations depend on this)
 */
expression::expression(unsigned int n,              // n. inputs
                   unsigned int m,                  // n. outputs
                   unsigned int r,                  // n. rows
                   unsigned int c,                  // n. columns
                   unsigned int l,                  // n. levels-back
                   function_set f,                  // functions
                   unsigned int seed                // seed for the pseudo-random numbers
                   ) : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_f(f), m_lb((3 * m_r * m_c) + m_m, 0), m_ub((3 * m_r * m_c) + m_m, 0), m_x((3 * m_r * m_c) + m_m, 0), m_e(seed)
{

    if (n == 0) throw input_error("Number of inputs is 0");
    if (m == 0) throw input_error("Number of outputs is 0");
    if (c == 0) throw input_error("Number of columns is 0");
    if (r == 0) throw input_error("Number of rows is 0");
    if (l == 0) throw input_error("Number of level-backs is 0");
    if (f.size()==0) throw input_error("Number of basis functions is 0");

    // Bounds for the function genes
    for (auto i = 0u; i < (3 * m_r * m_c); i+=3) {
        m_ub[i] = f.size() - 1;
    }

    // Bounds for the output genes
    for (auto i = 3u * m_r * m_c; i < m_ub.size(); ++i) {
        m_ub[i] = m_n + m_r * m_c - 1;
        if (m_l <= m_c) {
            m_lb[i] = m_n + m_r * (m_c - m_l);
        }
    }

    // Bounds for the node connection genes 
    for (auto i = 0u; i < m_c; ++i) {
        for (auto j = 0u; j < m_r; ++j) {
            m_ub[((i * m_r) + j) * 3 + 1] = m_n + i * m_r - 1;
            m_ub[((i * m_r) + j) * 3 + 2] = m_n + i * m_r - 1;
            if (i >= m_l) {
                m_lb[((i * m_r) + j) * 3 + 1] = m_n + m_r * (i - m_l);
                m_lb[((i * m_r) + j) * 3 + 2] = m_n + m_r * (i - m_l);
            }
        }
    }

    // We generate a random expression
    for (auto i = 0u; i < m_x.size(); ++i)
    {
        m_x[i] = std::uniform_int_distribution<unsigned int>(m_lb[i], m_ub[i])(m_e);
    }
    update_active();
}

/// Sets the chromosome
/** Sets a new chromosome as genotype for the expression and updates the active nodes and active genes information
 *
 * \param[in] x The new cromosome
 *
 * @throw dcgp::input_error if the chromosome is incompatible with the expression (n.inputs, n.outputs, levels-back, etc.)
 */
void expression::set(const std::vector<unsigned int>& x)
{
    if(!is_valid(x))
    {
        throw input_error("Chromosome is incompatible");
    }
    m_x = x;
    update_active();
}


inline unsigned int factorial(unsigned int n)
{
    if (n==0) return 1;
    unsigned int ret = 1;
    for(auto i = 1u; i <= n; ++i)
        ret *= i;
    return ret;
}

/// Accesses the derivatives of the expression
/** 
 * This method extracts a given derivative from a DA expression.
 *
 * \param[in] wrt array of variable indeces (0,1 ..., m_n) with respect to which to differentiate
 * \param[in] exp std::vector containing the DA expansion around the point where we want to compute the derivative
 *
 * @returns std::vector<double> containing the requested derivative of exp at the expansion point
 *
 * @throw dcgp::input_error 
 */
std::vector<double> expression::differentiate(const std::vector<unsigned int>& wrt, std::vector<DACE::DA> exp) const
{
    for (auto i : wrt)
    {
        if (i>=m_n)
        {
            throw input_error("Derivative id is larger than the independent variable number");
        }
        for (auto j = 0u; j<m_m; ++j)
        {
            exp[j] = exp[j].deriv(i+1u);
        }
    }

    std::vector<double> res(m_m);
    for (auto j = 0u; j<m_m; ++j)
    {
        res[j] = cons(exp[j]);
    }

    return res;
}

/// Computes the derivatives of the expression
/** 
 * Using differential algebra this method returns a DA expression containing all derivative at a given point.
 *
 * \param[in] wrt array of variable indeces (0,1 ..., m_n) with respect to which to differentiate
 * \param[in] in std::vector containing the point coordinates we want the derivative to be computed at
 *
 * @returns std::vector<double> containing the requested derivative of f at the point in
 *
 * @throw dcgp::input_error 
 */
std::vector<DACE::DA> expression::differentiate(const std::vector<double>& in) const
{  
    if (in.size() != m_n)
    {
        throw input_error("Input size is incompatible");
    }

    std::vector<DACE::DA> inDA(m_n);
    for (auto i = 0u; i < m_n; ++i)
    {
        inDA[i] = in[i] + DACE::DA(i+1);
    }
    return this->operator()(inDA);
}

/// Computes the derivatives of the expression
/** 
 * Using differential algebra this method returns a given derivative at a given point.
 * If several derivatives are needed at the same point, cache the result of differentiate(in)
 * and use differentiate(wrt,exp) on it repeatedly with different values for wrt.
 *
 * \param[in] wrt array of variable indeces (0,1 ..., m_n) with respect to which to differentiate
 * \param[in] in std::vector containing the point coordinates we want the derivative to be computed at
 *
 * @returns std::vector<double> containing the requested derivative of f at the point in
 *
 * @throw dcgp::input_error 
 */
std::vector<double> expression::differentiate(const std::vector<unsigned int>& wrt, const std::vector<double>& in) const
{  
    return differentiate(wrt,differentiate(in));
}

/// Mutates one of the active genes
/** 
 * Mutates exactly one of the active genes
 */
void expression::mutate_active()
{
    unsigned int idx = std::uniform_int_distribution<unsigned int>(0, m_active_genes.size() - 1)(m_e);
    idx = m_active_genes[idx];

    if (m_lb[idx]<m_ub[idx]) // if only one value is allowed for the gene, then we will not do anything as mutation does not apply
    {
        unsigned int new_value = UINT_MAX;
        do 
        {
            new_value = std::uniform_int_distribution<unsigned int>(m_lb[idx], m_ub[idx])(m_e);
        } while (new_value == m_x[idx]);
        m_x[idx] = new_value;
        update_active();
    }
}

/// Validity of a chromosome
/** 
 * Checks if a chromosome (i.e. a sequence of integers) is a valid expression
 * by checking its length and the bounds
 *
 * \param[in] x chromosome 
 */
bool expression::is_valid(const std::vector<unsigned int>& x) const
{
    // Checking for length
    if (x.size() != m_lb.size()) {
        return false;
    }

    // Checking for bounds on all cenes
    for (auto i = 0u; i < x.size(); ++i) {
        if ((x[i] > m_ub[i]) || (x[i] < m_lb[i])) {
            return false;
        }
    }
    return true;
}


/// Updates the m_active_genes and m_active_nodes data member
void expression::update_active()
{
    assert(m_x.size() == m_lb.size());

    // First we update the active nodes
    std::vector<unsigned int> current(m_m), next;
    m_active_nodes.clear();
    // At the beginning current contains only the output nodes connections
    for (auto i = 0u; i < m_m; ++i) {
        current[i] = m_x[3 * m_r * m_c + i];
    }
    do
    {
        m_active_nodes.insert(m_active_nodes.end(), current.begin(), current.end());

        for (auto node_id : current)
        {
            if (node_id >=m_n) // we insert the input nodes connections as they do not have any
            {
                unsigned int idx = (node_id - m_n) * 3;
                switch( m_f[m_x[idx]].m_type )
                {
                    case basis_function::FN_BINARY:
                        next.push_back(m_x[idx+2]);
                        // fall through
                    case basis_function::FN_UNARY:
                        next.push_back(m_x[idx+1]);
                        // fall through
                    case basis_function::FN_CONST:
                        break;
                }
            }
            else{
                m_active_nodes.push_back(node_id);
            }
        }
        // We remove duplicates to avoid processing them and thus having a 2^N complexity
        std::sort( next.begin(), next.end() );
        next.erase( std::unique( next.begin(), next.end() ), next.end() );
        current = next;
        next.clear();
    } while (current.size() > 0);

    // We remove duplicates and keep m_active_nodes sorted
    std::sort( m_active_nodes.begin(), m_active_nodes.end() );
    m_active_nodes.erase( std::unique( m_active_nodes.begin(), m_active_nodes.end() ), m_active_nodes.end() );

    // Then the active genes
    m_active_genes.clear();
    for (auto i = 0u; i<m_active_nodes.size(); ++i) 
    {
        if (m_active_nodes[i] >= m_n) 
        {
            unsigned int idx = (m_active_nodes[i] - m_n) * 3;
            m_active_genes.push_back(idx);
            m_active_genes.push_back(idx + 1);
            m_active_genes.push_back(idx + 2);
        }
    }
    for (auto i = 0u; i<m_m; ++i) 
    {
        m_active_genes.push_back(m_r * m_c * 3 + i);
    }
}

/// Return human readable representation of the problem.
/**
 * Will return a formatted string containing a human readable representation of the class
 *
 * @return std::string containing a human-readable representation of the problem.
 */
std::string expression::human_readable() const
{
    std::ostringstream s;
    s << "d-CGP Expression:\n";
    s << "\tNumber of inputs:\t\t" << m_n << '\n';
    s << "\tNumber of outputs:\t\t" << m_m << '\n';
    s << "\tNumber of rows:\t\t\t" << m_r << '\n';
    s << "\tNumber of columns:\t\t" << m_c << '\n';
    s << "\tNumber of levels-back allowed:\t" << m_l << '\n';
    s << "\n\tResulting lower bounds:\t" << m_lb;
    s << "\n\tResulting upper bounds:\t" << m_ub << '\n';
    s << "\n\tCurrent expression (encoded):\t" << m_x << '\n';
    s << "\tActive nodes:\t\t\t" << m_active_nodes << '\n';
    s << "\tActive genes:\t\t\t" << m_active_genes << '\n';
    s << "\n\tFunction set:\t\t\t" << m_f << '\n';
    return s.str();
}

/// Overload stream operator for dcgp::expression
/**
 * Equivalent to printing expression::human_readable() to stream.
 *
 * @param[out] s std::ostream to which the problem will be streamed.
 * @param[in] p dcgp::expression to be inserted into the stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const expression &en)
{
    s << en.human_readable();
    return s;
}

} // end of namespace dcgp

