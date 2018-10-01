/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_EIGEN_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#if defined EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "custom_utilities/ublas_wrapper.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"


namespace Kratos
{

template <typename Scalar>
struct SpaceType;

template <>
struct SpaceType<double>
{
	using Global = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
	using Local = UblasSpace<double, Matrix, Vector>;
};

template <>
struct SpaceType<std::complex<double>>
{
	using Global = UblasSpace<std::complex<double>, ComplexCompressedMatrix, boost::numeric::ublas::vector<std::complex<double>>>;
	using Local = UblasSpace<std::complex<double>, ComplexMatrix, ComplexVector>;
};

template <typename scalar_t>
struct SolverType
{
    using TScalar = scalar_t;

    using TSparseMatrix = Eigen::SparseMatrix<scalar_t, Eigen::RowMajor, int>;

    using TGlobalSpace = typename SpaceType<scalar_t>::Global;

    using TLocalSpace = typename SpaceType<scalar_t>::Local;
};

template <typename scalar_t>
struct SparseLU : public SolverType<scalar_t>
{
    using TSolver = Eigen::SparseLU<typename SolverType<scalar_t>::TSparseMatrix>;

    static constexpr auto Name = "SparseLU";
};

template <typename scalar_t>
struct SparseQR : SolverType<scalar_t>
{
    using TSolver = Eigen::SparseQR<typename SolverType<scalar_t>::TSparseMatrix, Eigen::COLAMDOrdering<int>>;

    static constexpr auto Name = "SparseQR";
};

#if defined EIGEN_USE_MKL_ALL
template <typename scalar_t>
struct PardisoLLT : SolverType<scalar_t>
{
    using TSolver = Eigen::PardisoLLT<typename SolverType<scalar_t>::TSparseMatrix>;

    static constexpr auto Name = "PardisoLLT";
};

template <typename scalar_t>
struct PardisoLDLT : SolverType<scalar_t>
{
    using TSolver = Eigen::PardisoLDLT<typename SolverType<scalar_t>::TSparseMatrix>;

    static constexpr auto Name = "PardisoLDLT";
};

template <typename scalar_t>
struct PardisoLU : SolverType<scalar_t>
{
    using TSolver = Eigen::PardisoLU<typename SolverType<scalar_t>::TSparseMatrix>;

    static constexpr auto Name = "PardisoLU";
};
#endif

template <
    class TSolverType,
    class TSparseSpaceType = typename TSolverType::TGlobalSpace,
    class TDenseSpaceType = typename TSolverType::TLocalSpace,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class EigenDirectSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
    typename TSolverType::TSolver m_solver;

    EigenDirectSolver &operator=(const EigenDirectSolver &Other);

    EigenDirectSolver(const EigenDirectSolver &Other);

  public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenDirectSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigenDirectSolver() {}

    EigenDirectSolver(Parameters settings) : BaseType(settings) {}

    ~EigenDirectSolver() override {}

    /**
     * This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        UblasWrapper<DataType> a_wrapper(rA);

        const auto& a = a_wrapper.matrix();

        m_solver.compute(a);

        KRATOS_ERROR_IF(m_solver.info() != Eigen::Success) << "Decomposition failed!" << std::endl;
    }

    /**
     * This function actually performs the solution work, eventually taking advantage of what was done before in the
     * Initialize and InitializeSolutionStep functions.
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1> > x(std::begin(rX.data()), rX.size());
        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1> > b(std::begin(rB.data()), rB.size());

        x = m_solver.solve(b);

        KRATOS_ERROR_IF(m_solver.info() != Eigen::Success) << "Solving failed!" << std::endl;
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    bool Solve(SparseMatrixType &rA, VectorType &rX, VectorType &rB) override
    {
        InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);

        return true;
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution matrix
     * @param rB Right hand side matrix
     * @return true if solution found, otherwise false
     */
    bool Solve(SparseMatrixType &rA, DenseMatrixType &rX, DenseMatrixType &rB) override
    {
        return false;
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "EigenDirectSolver<" << TSolverType::Name << "> finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }
}; // class EigenDirectSolver

/**
 * input stream function
 */
template<
    class TSolverType,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType>
inline std::istream &operator>>(
    std::istream &rIStream,
    EigenDirectSolver<TSolverType, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<
    class TSolverType,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType
    >
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EigenDirectSolver<TSolverType, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_SOLVER_H_INCLUDED)
