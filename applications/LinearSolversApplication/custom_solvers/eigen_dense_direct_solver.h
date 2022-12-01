/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_EIGEN_DENSE_DIRECT_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_DIRECT_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "custom_factories/dense_linear_solver_factory.h"
#include "custom_solvers/eigen_direct_solver.h"

namespace Kratos {

template <
    class TSolverType,
    class TSparseSpaceType = typename SpaceType<typename TSolverType::Scalar>::Local,
    class TDenseSpaceType = typename SpaceType<typename TSolverType::Scalar>::Local,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class EigenDenseDirectSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
private:
    TSolverType m_solver;

    EigenDenseDirectSolver &operator=(const EigenDenseDirectSolver &Other);

    EigenDenseDirectSolver(const EigenDenseDirectSolver &Other);

public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenDenseDirectSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType MatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigenDenseDirectSolver() {}

    EigenDenseDirectSolver(Parameters settings) : BaseType(settings)
    {
        m_solver.Initialize(settings);
    }

    ~EigenDenseDirectSolver() override {}

    /**
     * This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void InitializeSolutionStep(MatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Eigen::Map<Kratos::EigenDynamicMatrix<DataType>> A(rA.data().begin(), rA.size1(), rA.size2());

        const bool success = m_solver.Compute(A);

        KRATOS_ERROR_IF(!success) << "Decomposition failed!" << std::endl;
    }

    /**
     * This function actually performs the solution work, eventually taking advantage of what was done before in the
     * Initialize and InitializeSolutionStep functions.
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void PerformSolutionStep(MatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Eigen::Map<Kratos::EigenDynamicVector<DataType>> x(rX.data().begin(), rX.size());
        Eigen::Map<Kratos::EigenDynamicVector<DataType>> b(rB.data().begin(), rB.size());

        const bool success = m_solver.Solve(b, x);

        KRATOS_ERROR_IF(!success) << "Solving failed!\n" << m_solver.GetSolverErrorMessages() << std::endl;
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    bool Solve(MatrixType &rA, VectorType &rX, VectorType &rB) override
    {
        InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);

        return true;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system AX=B and puts the result on SystemMatrix& rX.
    @param rA. System matrix
    @param rX. Solution matrix.
     @param rB. Right hand side matrix.
    */
    bool Solve(MatrixType& rA, MatrixType& rX, MatrixType& rB) override
    {
        VectorType dummy;
        InitializeSolutionStep(rA, dummy, dummy);

        Eigen::Map<Kratos::EigenDynamicMatrix<DataType>> X(rX.data().begin(), rX.size1(), rX.size2());
        Eigen::Map<Kratos::EigenDynamicMatrix<DataType>> B(rB.data().begin(), rB.size1(), rB.size2());

        const bool success = m_solver.SolveMultiple(B, X);

        KRATOS_ERROR_IF(!success) << "Solving failed!\n" << m_solver.GetSolverErrorMessages() << std::endl;

        return true;
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

    static DenseLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenDenseDirectSolver> Factory()
    {
        return DenseLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenDenseDirectSolver>();
    }
}; // class EigenDenseDirectSolver

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
    EigenDenseDirectSolver<TSolverType, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
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
    const EigenDenseDirectSolver<TSolverType, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_DENSE_DIRECT_SOLVER_H_INCLUDED)
