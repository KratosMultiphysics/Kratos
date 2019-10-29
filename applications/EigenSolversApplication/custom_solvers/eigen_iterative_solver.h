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

// Project includes
#include "includes/define.h"
#include "custom_utilities/ublas_wrapper.h"
#include "factories/standard_linear_solver_factory.h"
#include "includes/ublas_complex_interface.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/iterative_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos {

template <typename Scalar>
struct SpaceType;

template <>
struct SpaceType<double>
{
    using Global = UblasSpace<double, CompressedMatrix, Vector>;
    using Local = UblasSpace<double, Matrix, Vector>;
};

template <>
struct SpaceType<std::complex<double>>
{
    using Global = UblasSpace<std::complex<double>, ComplexCompressedMatrix, ComplexVector>;
    using Local = UblasSpace<std::complex<double>, ComplexMatrix, ComplexVector>;
};

template <
    class TSolverType,
    class TSparseSpaceType = typename SpaceType<TSolverType::Scalar>::Global,
    class TDenseSpaceType = typename SpaceType<TSolverType::Scalar>::Local,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class EigenIterativeSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
private:
    typename TSolverType m_solver;

    EigenIterativeSolver &operator=(const EigenIterativeSolver &Other);

    EigenIterativeSolver(const EigenIterativeSolver &Other);

public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenIterativeSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef StandardLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenIterativeSolver> FactoryType;

    EigenIterativeSolver() {}

    EigenIterativeSolver(Parameters settings) : BaseType(settings)
    {
        m_solver.Initialize(settings);
    }

    ~EigenIterativeSolver() override {}

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

        const bool success = m_solver.Compute(a);

        KRATOS_ERROR_IF(!success) << "Decomposition failed!" << std::endl;
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
        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1>> x(rX.data().begin(), rX.size());
        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1>> b(rB.data().begin(), rB.size());

        const bool success = m_solver.Solve(b, x);

        KRATOS_ERROR_IF(!success) << "Decomposition failed!" << std::endl;
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
        m_solver.PrintInfo(rOStream);
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

    static StandardLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenIterativeSolver> Factory()
    {
        return StandardLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenIterativeSolver>();
    }
}; // class EigenIterativeSolver

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
    EigenIterativeSolver<TSolverType, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
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
    const EigenIterativeSolver<TSolverType, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_SOLVER_H_INCLUDED)
