//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//

#include <array>
#include <iomanip>
#include <cmath>
#include <ostream>
#include <iostream>

/**
 * This namespace contains classes and function mimicking the necessary
 * behaviour of Kratos in order to run the generated code without having
 * to modify it.
 */
namespace TestCompressibleNavierStokesSymbolic {

/** Matrix-like class with only:
 * - storage manipulation access
 * - operator<<
 * - almost equality check
 */
template<std::size_t TRows, std::size_t TCols>
class Matrix
{
public:
    using data_t = double;
    using index_t = std::size_t;

    Matrix(data_t fill = 0)
    {
        this->fill(fill);
    };

    Matrix& fill(data_t fill)
    {
        for(auto & row: m_data)
        {
            std::fill(row.begin(), row.end(), fill);
        }  
        return *this;
    }

    /**
     * Performs three cheks on a matrix:
     * - Ensures that no entries are NaN
     * - Ensures that no entries are inf
     * - Ensures almost-equality with a reference matrix
     *
     * @param[in]  reference  The reference matrix to compare to
     * @param[in]  rel_error  The maximum allowable relative error
     * 
     * @return An error code. It is 0 if no errors, 1 if there are
     * The errors are also be printed to std::err
     */
    template<index_t R, index_t C>
    int validate(Matrix<R,C> const& reference, const double rel_error = 1e-6) const
    {
        static_assert(R==TRows, "Diferent number of rows!");
        static_assert(C==TCols, "Diferent number of cols!");

        int result = 0;
        for(index_t i=0; i < R*C; ++i)
        {
            if(std::isnan((*this)[i]))
            {
                std::cerr << "\nEntry #" << i << " is NaN";
                result = 1;
            }
            else if(std::isinf((*this)[i]))
            {
                std::cerr << "\nEntry #" << i << " is is infinite";
                result = 1;   
            }
        }

        for(index_t i=0; i < R*C; ++i)
        {
            const data_t value = (*this)[i];
            const data_t ref = reference[i];

            if(value==0 && ref==0) continue;

            const data_t diff = value - ref;
            const data_t max = std::max(std::abs(value), std::abs(ref));

            const data_t delta = rel_error * max;
            const int precision = 1 - log10(rel_error);
    
            if(std::abs(diff) > delta)
            {
                std::cerr << std::setprecision(precision)
                          << "\nMismatch in entry #" << i <<":\t" 
                          << (*this)[i] << " != " << reference[i] 
                          << "\t(delta = " 
                          << std::setprecision(2) << delta
                          << ").\tDiff = " << diff
                          << "\tRel diff = " << diff/max;
                result = 1;
            }
        }

        return result;
    }

    constexpr data_t const& operator[](const index_t i) const
    {
        return (*this)(i);
    }

    data_t & operator[](const index_t i)
    {
        return (*this)(i);
    }

    constexpr data_t const& operator()(const index_t i, const index_t j) const
    {
        return m_data[i][j];
    }

    constexpr data_t const& operator()(const index_t i) const
    {
        return m_data[i / TCols][i % TCols];
    }

    data_t & operator()(const index_t i, const index_t j)
    {
        return m_data[i][j];
    }

    data_t & operator()(const index_t i)
    {
        return m_data[i / TCols][i % TCols];
    }

    friend std::ostream& operator<<(std::ostream& os, Matrix const& m)
    {
        for(auto const& row: m.m_data)
        {
            os << "[  " << row[0];
            for(index_t i=1; i<TCols; ++i)
            {
                os << ", " << row[i];
            }
            os << "    ]\n";
        }
        os << "\n";

        return os;
    }

private:
    std::array<std::array<data_t, TCols>, TRows> m_data;
};


template<std::size_t TRows>
using Vector = Matrix<TRows, 1>;

/**
 * @brief Shape functions in isoparametric shapes
 */
template<std::size_t TDim, std::size_t TNodes>
void ShapeFunctions(Vector<TNodes> & N, Matrix<TNodes, TDim> & DN_DX, double x, double y);

template<>
inline void ShapeFunctions(Vector<4> & N, Matrix<4, 2> & DN_DX, double x, double y)
{
    N[0] = (1-x)*(1-y)/4;
    N[1] = (1+x)*(1-y)/4;
    N[2] = (1+x)*(1+y)/4;
    N[3] = (1-x)*(1+y)/4;

    DN_DX(0,0) = -(1-y)/4;      DN_DX(0,1) = -(1-x)/4;
    DN_DX(1,0) =  (1-y)/4;      DN_DX(1,1) = -(1+x)/4;
    DN_DX(2,0) =  (1+y)/4;      DN_DX(2,1) =  (1+x)/4;
    DN_DX(3,0) = -(1+y)/4;      DN_DX(3,1) =  (1-x)/4;
}

template<>
inline void ShapeFunctions(Vector<3> & N, Matrix<3, 2> & DN_DX, double x, double y)
{
    N[0] = -(x+y)/2;
    N[1] = (1+x)/2;;
    N[2] = (1+y)/2;

    DN_DX(0,0) = -0.5;      DN_DX(0,1) = -0.5;
    DN_DX(1,0) =  0.5;      DN_DX(1,1) =  0;
    DN_DX(2,0) =  0;        DN_DX(2,1) =  0.5;
}

/**
 * Struct containing data mimicking
 * CompressibleNavierStokesExplicit::ElementDataStruct.
 * 
 * By default filled with Rankine-Hugoniot jump condition.
 * Difusive magnitudes are made-up but within the correct order of magnitude.
 */
template<unsigned int Tdim, unsigned int Tnodes, unsigned int Tblocksize = Tdim+2>
struct ElementDataT
{
    Matrix<Tblocksize, Tnodes> U;
    Matrix<Tblocksize, Tnodes> dUdt;
    
    Vector<Tnodes> m_ext;
    Vector<Tnodes> r_ext;
    Matrix<Tnodes, Tdim> f_ext;

    Vector<Tnodes> alpha_sc_nodes;
    Vector<Tnodes> beta_sc_nodes;
    Vector<Tnodes> lamb_sc_nodes;
    Vector<Tnodes> mu_sc_nodes;
    double alpha, beta, lambda, mu;

    Matrix<Tnodes, Tblocksize> ResProj;

    double h;
    double gamma, c_v;

    ElementDataT();
};

/// Quad
template<>
inline ElementDataT<2, 4>::ElementDataT()
{
    constexpr double rho_0 = 1.16927;
    constexpr double rho_1 = 1.46426;

    constexpr double mom = 467.707;

    constexpr double et_0 = 346854;
    constexpr double et_1 = 422234;

    this->U(0, 0) = rho_0;
    this->U(1, 0) = rho_1;
    this->U(2, 0) = rho_1;
    this->U(3, 0) = rho_0;

    this->U(0, 1) = mom;
    this->U(1, 1) = mom;
    this->U(2, 1) = mom;
    this->U(3, 1) = mom;

    this->U(0, 3) = et_0;
    this->U(1, 3) = et_1;
    this->U(2, 3) = et_1;
    this->U(3, 3) = et_0;

    this->alpha_sc_nodes.fill(1.5e-4);
    this->beta_sc_nodes.fill(2.8e-5);
    this->lamb_sc_nodes.fill(1.3e-7);
    this->mu_sc_nodes.fill(2.3e-6);

    this->alpha = 0;
    this->beta = 1.13e-4;
    this->lambda = 6.84e-6;
    this->mu = 1.26e-4;

    this->gamma = 1.4;
    this->c_v = 722.14;
    this->h = 2.0;
}

/// Triangle
template<>
inline ElementDataT<2, 3>::ElementDataT()
{
    constexpr double rho_0 = 1.16927;
    constexpr double rho_1 = 1.46426;

    constexpr double mom = 467.707;

    constexpr double et_0 = 346854;
    constexpr double et_1 = 422234;

    this->U(0, 0) = rho_0;
    this->U(1, 0) = rho_1;
    this->U(2, 0) = rho_1;

    this->U(0, 1) = mom;
    this->U(1, 1) = mom;
    this->U(2, 1) = mom;

    this->U(0, 3) = et_0;
    this->U(1, 3) = et_1;
    this->U(2, 3) = et_1;

    this->alpha_sc_nodes.fill(1.5e-4);
    this->beta_sc_nodes.fill(2.8e-5);
    this->lamb_sc_nodes.fill(1.3e-7);
    this->mu_sc_nodes.fill(2.3e-6);

    this->alpha = 0;
    this->beta = 1.13e-4;
    this->lambda = 6.84e-6;
    this->mu = 1.26e-4;

    this->gamma = 1.4;
    this->c_v = 722.14;
    this->h = 2.0;
}

}
