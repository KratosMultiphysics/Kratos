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
#include "boost/smart_ptr.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#if defined EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/direct_solver.h"

#include <chrono>
using namespace std::chrono;

namespace Kratos
{

struct SolverType
{
    using TSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
};

struct SparseLU : SolverType
{
    using TSolver = Eigen::SparseLU<TSparseMatrix>;

    static constexpr auto Name = "SparseLU";
};

#if defined EIGEN_USE_MKL_ALL
struct PardisoLLT : SolverType
{
    using TSolver = Eigen::PardisoLLT<TSparseMatrix>;

    static constexpr auto Name = "PardisoLLT";
};

struct PardisoLDLT : SolverType
{
    using TSolver = Eigen::PardisoLDLT<TSparseMatrix>;

    static constexpr auto Name = "PardisoLDLT";
};

struct PardisoLU : SolverType
{
    using TSolver = Eigen::PardisoLU<TSparseMatrix>;

    static constexpr auto Name = "PardisoLU";
};
#endif

template <
    class TSolver,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class EigenDirectSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
    EigenDirectSolver &operator=(const EigenDirectSolver &Other);

    EigenDirectSolver(const EigenDirectSolver &Other);

  public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenDirectSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigenDirectSolver() {}

    EigenDirectSolver(Parameters settings) : BaseType(settings) {}

    ~EigenDirectSolver() override {}

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    bool Solve(SparseMatrixType &rA, VectorType &rX, VectorType &rB) override
    {
        std::vector<int> index1_vector(rA.index1_data().size());
        std::vector<int> index2_vector(rA.index2_data().size());

        for (size_t i = 0; i < rA.index1_data().size(); i++) {
            index1_vector[i] = (int)rA.index1_data()[i];
        }

        for (size_t i = 0; i < rA.index2_data().size(); i++) {
            index2_vector[i] = (int)rA.index2_data()[i];
        }

        Eigen::Map<typename TSolver::TSparseMatrix> a(rA.size1(), rA.size2(), rA.nnz(), index1_vector.data(), index2_vector.data(), rA.value_data().begin());
        Eigen::Map<Eigen::VectorXd> x(rX.data().begin(), rX.size());
        Eigen::Map<Eigen::VectorXd> b(rB.data().begin(), rB.size());

        typename TSolver::TSolver solver;
        solver.compute(a);
        x = solver.solve(b);

        bool success = (solver.info() == Eigen::Success);

        return success;
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
        rOStream << "EigenDirectSolver<" << TSolver::Name << "> finished.";
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
    class TSolver,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType>
inline std::istream &operator>>(
    std::istream &rIStream,
    EigenDirectSolver<TSolver, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<
    class TSolver,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType
    >
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EigenDirectSolver<TSolver, TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_SOLVER_H_INCLUDED)