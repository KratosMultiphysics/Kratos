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

// System includes
#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>
#include <sstream>

// External includes

// Project includes

///@name Globals
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

enum class TestResult { SUCCESS, FAILURE };

///@}
///@name  Functions
///@{

/**
 * Syntactic sugar to accumulate test results
 */
constexpr TestResult operator+=(TestResult & lhs, const TestResult& rhs)
{
    return lhs = (lhs == TestResult::SUCCESS && rhs == TestResult::SUCCESS)
        ? TestResult::SUCCESS
        : TestResult::FAILURE;
}

///@}
///@name Classes
///@{


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
     * Validates that all entries are real valued (i.e. neither inf nor nan)
     */
    bool ValuesAreReal(std::ostream& error_stream = std::cerr) const
    {
        bool result = true;
        for(IndexType i = 0; i < TRows*TCols; ++i)
        {
            if(std::isnan((*this)[i]))
            {
                error_stream << "\nEntry #" << i << " is NaN";
                result = false;
            }
            else if(std::isinf((*this)[i]))
            {
                error_stream << "\nEntry #" << i << " is is infinite";
                result = false;
            }
        }
        return result;
    }

    /**
     * Returns whether the matrix is almost equal to the reference,
     * within relative tolerance.
     */
    bool AlmostEqual(
        const Matrix& reference,
        const DataType rel_tolerance,
        std::ostream& error_stream = std::cerr) const
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
                error_stream << std::setprecision(precision)
                             << "Mismatch in entry #"
                             << std::setw(3) << i
                             <<":"
                             << std::setw(10) << (*this)[i]
                             << " != "
                             << std::setw(10) << reference[i]
                             << " | "
                             << "Diff = "
                             << std::setw(10) << diff
                             << " | "
                             << "Rel diff = " << diff/max << '\n';
                result = false;
            }
        }

        return result;
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
void ShapeFunctions(Vector<TNodes>& N, Matrix<TNodes, TDim>& DN_DX, double x, double y = 0, double z = 0);

/**
 * Triangle shape functions
 */
template<>
inline void ShapeFunctions(Vector<3>& N, Matrix<3, 2>& DN_DX, double x, double y, double)
{
    N[0] = -(x+y)/2;
    N[1] = (1+x)/2;
    N[2] = (1+y)/2;

    DN_DX(0,0) = -0.5;      DN_DX(0,1) = -0.5;
    DN_DX(1,0) =  0.5;      DN_DX(1,1) =  0;
    DN_DX(2,0) =  0;        DN_DX(2,1) =  0.5;
}

/**
 * Quadrilateral shape functions
 */
template<>
inline void ShapeFunctions(Vector<4>& N, Matrix<4, 2>& DN_DX, double x, double y, double)
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
 * Tetrahedron shape functions
 */
template<>
inline void ShapeFunctions(Vector<4>& N, Matrix<4, 3>& DN_DX, double x, double y, double z)
{
    N[0] = -(1+x+y+z)/3;
    N[1] = (1+x)/2;
    N[2] = (1+y)/2;
    N[3] = (1+z)/2;

    DN_DX(0,0) = -1/3.0;    DN_DX(0,1) = -1/3.0;    DN_DX(0,2) = -1/3.0;
    DN_DX(1,0) =  0.5;      DN_DX(1,1) =  0;        DN_DX(1,2) = 0;
    DN_DX(2,0) =  0;        DN_DX(2,1) =  0.5;      DN_DX(2,2) = 0;
    DN_DX(3,0) =  0;        DN_DX(3,1) =  0;        DN_DX(3,2) = 0.5;
}

/**
 * Struct containing data mimicking
 * CompressibleNavierStokesExplicit::ElementDataStruct.
 *
 */
template<std::size_t Tdim, std::size_t Tnodes>
struct ElementData
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
};

/*********************************************************
 *                                                       *
 *                  Testing utilities                    *
 *                                                       *
 *********************************************************/

/**
 * Template to be used by CheckSubstitutionResult.
 *
 * Specializations must perform two tasks:
 *  - Ensures the state of `result` is valid
 *  - Ensures the state of `result` is within tolerance of `reference`
 *  - Any errors must be reported to error_stream
 */
template<typename T, typename U>
TestResult Validate(
    const T& result,
    const U& reference,
    const double rel_tolerance,
    std::ostream& error_stream = std::cerr);

/**
 * Specialization for matrices:
 *
 * Performs three cheks on a matrix:
 * - Ensures that no entries are NaN
 * - Ensures that no entries are inf
 * - Ensures almost-equality with a reference matrix
 *
 * @param  result         The matrix to test
 * @param  reference      The reference matrix to compare to
 * @param  rel_tolerance  The maximum allowable relative error
 * @param  error_stream   The stream to print the errors into
 *
 * @return The result of the test
 */
template<std::size_t TRows, std::size_t TCols>
inline TestResult Validate(
    const Matrix<TRows, TCols>& result,
    const Matrix<TRows, TCols>& reference,
    const double rel_tolerance,
    std::ostream& error_stream = std::cerr)
{
    if(!result.ValuesAreReal(error_stream))
    {
        return TestResult::FAILURE;
    }

    if(!result.AlmostEqual(reference, rel_tolerance, error_stream))
    {
        return TestResult::FAILURE;
    }

    return TestResult::SUCCESS;
}

/**
 * Boilerplate to check that the obtained result is equal to the reference,
 * with some information prints to stdout and error reports to stderr.
 *
 * A specialization Validate<T,U> must exist for the check to be performed
 */
template<typename T, typename U>
inline TestResult CheckSubstitutionResult(
    const std::string test_name,
    const T& result,
    const U& expected,
    const bool print_results)
{
    std::cout << "Testing result of //substitute_"<< test_name << ":";

    std::stringstream error_stream;

    const TestResult test_result = Validate(result, expected, 1e-6, error_stream);

    switch(test_result)
    {
        case TestResult::SUCCESS:
            std::cout << " Passed" << std::endl;
            break;
        case TestResult::FAILURE:
            std::cout << " Failed" << std::endl;
            std::cerr << "Errors in result of //substitute_"<< test_name << ":\n";
            std::cerr << error_stream.str() << std::endl;
    }

    if(print_results)
    {
        std::cout << "Results:\n" << std::setprecision(20) << result << std::endl;
    }

    return test_result;
}