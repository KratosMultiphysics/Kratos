#include <array>
#include <iomanip>
#include <cmath>
#include <ostream>
#include <iostream>

namespace TestCompressibleNavierStokesSymbolic {

template<std::size_t TDim>
struct Dofs {
    static constexpr std::size_t RHO = 0;
    static constexpr std::size_t MOM_X = 1;
    static constexpr std::size_t MOM_Y = 2;
    static constexpr std::size_t MOM_Z = 3;
    static constexpr std::size_t E_TOT = TDim + 1;
};

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

    template<index_t R, index_t C>
    int validate(Matrix<R,C> const& reference, const double delta = 1e-8) const
    {
        static_assert(R==TRows, "Diferent number of rows!");
        static_assert(C==TCols, "Diferent number of cols!");

        const int precision = 1 - log10(delta);

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
            if(std::abs((*this)[i] - reference[i]) > delta)
            {
                std::cerr << std::setprecision(precision)
                          << "\nMismatch in entry #" << i <<": " 
                          << (*this)[i] << " != " << reference[i] 
                          << " within delta<" 
                          << std::setprecision(2) << delta;
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


inline void QuadShapeFunctions(Vector<4> & N, Matrix<4, 2> & DN_DX, double x, double y)
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

template<unsigned int Tnodes, unsigned int Tdim, unsigned int Tblocksize>
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
};


ElementDataT<4, 2, 4> RankineHugoniotQuadData()
{
    ElementDataT<4,2,4> data;

    constexpr double rho_0 = 1.16927;
    constexpr double rho_1 = 1.46426;

    constexpr double mom = 467.707;

    constexpr double et_0 = 346854;
    constexpr double et_1 = 422234;

    using d = Dofs<2>;

    data.U(d::RHO, 0) = rho_0;
    data.U(d::RHO, 1) = rho_1;
    data.U(d::RHO, 2) = rho_1;
    data.U(d::RHO, 3) = rho_0;

    data.U(d::MOM_X, 0) = mom;
    data.U(d::MOM_X, 1) = mom;
    data.U(d::MOM_X, 2) = mom;
    data.U(d::MOM_X, 3) = mom;

    data.U(d::E_TOT, 0) = et_0;
    data.U(d::E_TOT, 1) = et_1;
    data.U(d::E_TOT, 2) = et_1;
    data.U(d::E_TOT, 3) = et_0;

    data.alpha_sc_nodes.fill(1.5e-4);
    data.beta_sc_nodes.fill(2.8e-5);
    data.lamb_sc_nodes.fill(1.3e-7);
    data.mu_sc_nodes.fill(2.3e-6);

    data.alpha = 0;
    data.beta = 1.13e-4;
    data.lambda = 6.84e-6;
    data.mu = 1.26e-4;

    data.gamma = 1.4;
    data.c_v = 722.14;
    data.h = 2.0;

    return data;
}

ElementDataT<4, 2, 4> SodQuadData()
{
    ElementDataT<4,2,4> data;

    constexpr double rho_0 = 1.0;
    constexpr double rho_1 = 0.125;

    constexpr double mom = 0.0;

    constexpr double et_0 = 2.5;
    constexpr double et_1 = 0.25;

    using d = Dofs<2>;

    data.U(d::RHO, 0) = rho_0;
    data.U(d::RHO, 1) = rho_1;
    data.U(d::RHO, 2) = rho_1;
    data.U(d::RHO, 3) = rho_0;

    data.U(d::MOM_X, 0) = mom;
    data.U(d::MOM_X, 1) = mom;
    data.U(d::MOM_X, 2) = mom;
    data.U(d::MOM_X, 3) = mom;

    data.U(d::E_TOT, 0) = et_0;
    data.U(d::E_TOT, 1) = et_1;
    data.U(d::E_TOT, 2) = et_1;
    data.U(d::E_TOT, 3) = et_0;

    data.alpha_sc_nodes.fill(1.5e-4);
    data.beta_sc_nodes.fill(2.8e-5);
    data.lamb_sc_nodes.fill(1.3e-7);
    data.mu_sc_nodes.fill(2.3e-6);

    data.alpha = 0;
    data.beta = 1.13e-4;
    data.lambda = 6.84e-6;
    data.mu = 1.26e-4;

    data.gamma = 1.4;
    data.c_v = 722.14;
    data.h = 2.0;

    return data;
}

}