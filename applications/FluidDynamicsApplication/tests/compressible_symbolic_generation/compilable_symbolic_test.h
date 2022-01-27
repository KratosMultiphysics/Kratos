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
//
// This file contains classes and functions with two roles:
// - Testing utilities
// - Mimicking the behaviour of Kratos in order to run the generated
//   code without having to modify it.

#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>


enum class TestResult { SUCCESS, FAILURE };

constexpr TestResult operator+=(TestResult & lhs, const TestResult& rhs)
{
    return (lhs == TestResult::SUCCESS && rhs == TestResult::SUCCESS)
        ? TestResult::SUCCESS 
        : TestResult::FAILURE;
}

/**
 * Boilerplate to check that the obtained result is equal to the reference,
 * with some information prints to stdout.
 */
template<typename T>
inline TestResult CheckSubstitutionResult(
    const std::string test_name,
    const T& result,
    const T& expected,
    const bool print_results)
{
    std::cout << "Testing result of //substitute_"<< test_name << ":";

    TestResult test_result = result.Validate(expected);

    if(test_result == TestResult::SUCCESS)
    {
        std::cout << " OK";
    }

    std::cout << std::endl;

    if(print_results)
    {
        std::cout << "Results:\n" << std::setprecision(20) << result << std::endl;
    }

    return test_result;
}


/** Matrix-like class with only:
 * - storage manipulation access
 * - operator<<
 * - almost equality check
 */
template<std::size_t TRows, std::size_t TCols>
class Matrix
{
public:
    static_assert(TRows != 0, "Number of rows must be greater than 0");
    static_assert(TCols != 0, "Number of colums must be greater than 0");

    ///@name Type Definitions
    ///@{
    
    using DataType = double;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Matrix(const DataType FillValue = 0)
    {
        this->Fill(FillValue);
    };

    /// Destructor.
    ~Matrix() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Fill(const DataType Value)
    {
        for(auto& row: mData)
        {
            std::fill(row.begin(), row.end(), Value);
        }
    }

    ///@}
    ///@name Access
    ///@{

    // Single index accessing

    constexpr DataType operator[](const IndexType Index) const
    {
        return mData[Index / TCols][Index % TCols];
    }

    DataType& operator[](const IndexType Index)
    {
        return mData[Index / TCols][Index % TCols];
    }

    constexpr DataType operator()(const IndexType Index) const
    {
        return (*this)[Index];
    }

    DataType& operator()(const IndexType Index)
    {
        return (*this)[Index];
    }

    // Multiple index accessing

    constexpr DataType operator()(const IndexType Row, const IndexType Col) const
    {
        return mData[Row][Col];
    }

    DataType& operator()(const IndexType Row, const IndexType Col)
    {
        return mData[Row][Col];
    }


    ///@}
    ///@name Inquiry
    ///@{

    /**
     * Performs three cheks on a matrix:
     * - Ensures that no entries are NaN
     * - Ensures that no entries are inf
     * - Ensures almost-equality with a reference matrix
     *
     * @param[in]  reference  The reference matrix to compare to
     * @param[in]  rel_tolerance  The maximum allowable relative error
     * 
     * @return An error code. It is 0 if no errors, 1 if there are
     * The errors are also be printed to std::err
     */
    TestResult Validate(const Matrix& reference, const double rel_tolerance = 1e-6) const
    {
        if(!ValuesAreReal()) return TestResult::FAILURE;
        if(!AlmostEqual(reference, rel_tolerance)) return TestResult::FAILURE;

        return TestResult::SUCCESS;
    }


    ///@}
    ///@name Input and output
    ///@{

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m)
    {
        for(const auto& r_row: m.mData)
        {
            os << "[  " << r_row[0];
            for(IndexType i=1; i<TCols; ++i)
            {
                os << ", " << r_row[i];
            }
            os << "    ]\n";
        }
        os << "\n";

        return os;
    }

    ///@}

private:
    ///@name Private static member variables
    ///@{


    ///@}
    ///@name Private member Variables
    ///@{

    std::array<std::array<DataType, TCols>, TRows> mData;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    /** 
     * Validates that all entries are real valued (i.e. neither inf nor nan)
     */
    bool ValuesAreReal() const
    {
        bool result = true;
        for(IndexType i = 0; i < TRows*TCols; ++i)
        {
            if(std::isnan((*this)[i]))
            {
                std::cerr << "\nEntry #" << i << " is NaN";
                result = false;
            }
            else if(std::isinf((*this)[i]))
            {
                std::cerr << "\nEntry #" << i << " is is infinite";
                result = false;   
            }
        }

        return result;
    }

    /** 
     * Returns whether the matrix is almost equal to the reference, 
     * within relative tolerance.
     */
    bool AlmostEqual(const Matrix& reference, const DataType rel_tolerance) const
    {
        bool result = true;

        for(IndexType i=0; i < TRows*TCols; ++i)
        {
            const DataType value = (*this)[i];
            const DataType ref = reference[i];

            if(value==0 && ref==0) continue;

            const DataType diff = value - ref;
            const DataType max = std::max(std::abs(value), std::abs(ref));

            const DataType delta = rel_tolerance * max;
            const int precision = 1 - std::log10(rel_tolerance);
    
            if(std::abs(diff) > delta)
            {
                std::cerr << std::setprecision(precision)
                          << "\nMismatch in entry #" << i <<":\t" 
                          << (*this)[i] << " != " << reference[i] 
                          << "\t(delta = " 
                          << std::setprecision(2) << delta
                          << ").\tDiff = " << diff
                          << "\tRel diff = " << diff/max;
                result = false;
            }
        }

        return result;
    }


    ///@}
    ///@name Private LifeCycle
    ///@{


    ///@}

};

///@}
///@name Type Definitions
///@{

template<std::size_t TRows>
using Vector = Matrix<TRows, 1>;

///@}
///@name Input and output
///@{


///@}


/**
 * @brief Shape functions in isoparametric space
 */
template<std::size_t TDim, std::size_t TNodes>
void ShapeFunctions(Vector<TNodes>& N, Matrix<TNodes, TDim>& DN_DX, double x, double y);

/**
 * Triangle shape functions
 */
template<>
inline void ShapeFunctions(Vector<3>& N, Matrix<3, 2>& DN_DX, double x, double y)
{
    N[0] = -(x+y)/2;
    N[1] = (1+x)/2;;
    N[2] = (1+y)/2;

    DN_DX(0,0) = -0.5;      DN_DX(0,1) = -0.5;
    DN_DX(1,0) =  0.5;      DN_DX(1,1) =  0;
    DN_DX(2,0) =  0;        DN_DX(2,1) =  0.5;
}

/**
 * Quadrilateral shape functions
 */
template<>
inline void ShapeFunctions(Vector<4>& N, Matrix<4, 2>& DN_DX, double x, double y)
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

/**
 * Struct containing data mimicking
 * CompressibleNavierStokesExplicit::ElementDataStruct.
 * 
 * By default filled with Rankine-Hugoniot jump condition.
 * Difusive magnitudes are made-up but within the correct order of magnitude.
 */
template<std::size_t Tdim, std::size_t Tnodes>
struct ElementDataT
{
    static constexpr std::size_t BlockSize = Tdim + 2;

    Matrix<BlockSize, Tnodes> U;
    Matrix<BlockSize, Tnodes> dUdt;
    
    Vector<Tnodes> m_ext;
    Vector<Tnodes> r_ext;
    Matrix<Tnodes, Tdim> f_ext;

    Vector<Tnodes> alpha_sc_nodes;
    Vector<Tnodes> beta_sc_nodes;
    Vector<Tnodes> lamb_sc_nodes;
    Vector<Tnodes> mu_sc_nodes;
    double alpha, beta, lambda, mu;

    Matrix<Tnodes, BlockSize> ResProj;

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

    this->alpha_sc_nodes.Fill(1.5e-4);
    this->beta_sc_nodes.Fill(2.8e-5);
    this->lamb_sc_nodes.Fill(1.3e-7);
    this->mu_sc_nodes.Fill(2.3e-6);

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

    this->alpha_sc_nodes.Fill(1.5e-4);
    this->beta_sc_nodes.Fill(2.8e-5);
    this->lamb_sc_nodes.Fill(1.3e-7);
    this->mu_sc_nodes.Fill(2.3e-6);

    this->alpha = 0;
    this->beta = 1.13e-4;
    this->lambda = 6.84e-6;
    this->mu = 1.26e-4;

    this->gamma = 1.4;
    this->c_v = 722.14;
    this->h = 2.0;
}