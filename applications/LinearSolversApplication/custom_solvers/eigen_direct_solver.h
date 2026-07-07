/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#pragma once

// System includes
#include <type_traits>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
#include "custom_utilities/ublas_wrapper.h"
#include "factories/standard_linear_solver_factory.h"
#include "includes/ublas_complex_interface.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "spaces/default_spaces.h"

namespace Kratos {

template <typename Scalar>
struct SpaceType;

// The real spaces follow the configure-time selected linear-algebra backend;
// the complex ones stay uBLAS (see spaces/default_spaces.h). UblasWrapper only
// relies on the compressed_matrix member surface, which both backends expose.
template <>
struct SpaceType<double>
{
    using Global = TDefaultSparseSpace<double>;
    using Local = TDefaultDenseSpace<double>;
};

template <>
struct SpaceType<std::complex<double>>
{
    using Global = TDefaultSparseSpace<std::complex<double>>;
    using Local = TDefaultDenseSpace<std::complex<double>>;
};

template <
    class TSolverType,
    class TSparseSpaceType = typename SpaceType<typename TSolverType::Scalar>::Global,
    class TDenseSpaceType = typename SpaceType<typename TSolverType::Scalar>::Local,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class EigenDirectSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
private:
    TSolverType m_solver;

    EigenDirectSolver &operator=(const EigenDirectSolver &Other);

    EigenDirectSolver(const EigenDirectSolver &Other);

    // The wrapper's map type matches what the solver's Compute expects, so the
    // copy fallback always produces the exact index type the solver needs.
    UblasWrapper<typename TSolverType::Scalar, typename TSolverType::SparseMatrix> m_a_wrapper;

public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenDirectSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef StandardLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenDirectSolver> FactoryType;

    EigenDirectSolver() {}

    EigenDirectSolver(Parameters settings) : BaseType(settings)
    {
        m_solver.Initialize(settings);
    }

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
        using SolverMatrixType = typename TSolverType::SparseMatrix;

        // Zero-copy fast path: when the system matrix is already an Eigen sparse
        // matrix whose storage index matches the solver's (the Eigen backend for
        // the pure-Eigen solvers), hand the solver a Map over its own CSR arrays
        // with no index rebuild. This is what makes an Eigen-backend solve cheaper
        // than uBLAS. Otherwise (uBLAS backend, or an index mismatch such as
        // MKL/Pardiso which requires int) copy once through UblasWrapper.
        constexpr bool is_native_eigen = requires (SparseMatrixType& rMatrix) {
            rMatrix.outerIndexPtr(); rMatrix.innerIndexPtr(); rMatrix.valuePtr();
            requires std::is_same_v<typename SparseMatrixType::StorageIndex,
                                    typename SolverMatrixType::StorageIndex>;
        };

        bool success;
        if constexpr (is_native_eigen) {
            KRATOS_DEBUG_ERROR_IF_NOT(rA.isCompressed())
                << "The Eigen system matrix must be in compressed (CSR) storage "
                   "before solving." << std::endl;
            const Eigen::Map<const SolverMatrixType> a(
                rA.rows(), rA.cols(), rA.nonZeros(),
                rA.outerIndexPtr(), rA.innerIndexPtr(), rA.valuePtr());
            success = m_solver.Compute(a);
        } else {
            m_a_wrapper = UblasWrapper<typename TSolverType::Scalar, SolverMatrixType>(rA);
            success = m_solver.Compute(m_a_wrapper.matrix());
        }

        KRATOS_ERROR_IF(!success) << "Decomposition failed!" << std::endl;
    }

    /**
     * This function actually performs the solution work, eventually taking advantage of what was done before in the
     * Initialize and InitializeSolutionStep functions.
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    bool PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        // uBLAS vectors expose the raw buffer through data().begin(), Eigen
        // ones directly through data()
        const auto get_buffer = [](VectorType& rVector) {
            if constexpr (requires { rVector.data().begin(); }) {
                return rVector.data().begin();
            } else {
                return rVector.data();
            }
        };
        Eigen::Map<Kratos::EigenDynamicVector<DataType>> x(get_buffer(rX), rX.size());
        Eigen::Map<Kratos::EigenDynamicVector<DataType>> b(get_buffer(rB), rB.size());

        const bool success = m_solver.Solve(b, x);

        KRATOS_ERROR_IF(!success) << "Solving failed!\n" << m_solver.GetSolverErrorMessages() << std::endl;
        return success;
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
        return PerformSolutionStep(rA, rX, rB);
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution matrix
     * @param rB Right hand side matrix
     * @return true if solution found, otherwise false
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        KRATOS_TRY

        const std::size_t system_size = TSparseSpaceType::Size1(rA);
        const std::size_t n_rhs = TDenseSpaceType::Size2(rB);
    
        typename TSparseSpaceType::VectorType solution, rhs;
        TSparseSpaceType::Resize(solution, TSparseSpaceType::Size1(rA));
        TSparseSpaceType::Resize(rhs, TSparseSpaceType::Size1(rA));
    
        KRATOS_ERROR_IF(TDenseSpaceType::Size1(rB) != system_size)
            << "Expecting the right hand side matrix to have "
            << system_size
            << " rows, but it has "
            << TDenseSpaceType::Size1(rB) << ".";
    
        TDenseSpaceType::Resize(rX, system_size, n_rhs);

        // Element-wise column extraction/insertion: works for any combination
        // of dense-matrix and system-vector backends and any scalar type
        // (unlike direct column() proxy assignment or the space GetColumn,
        // which is bound to the real dense matrix)
        const auto get_column = [](const std::size_t J, const DenseMatrixType& rM, VectorType& rColumn) {
            if (static_cast<std::size_t>(rColumn.size()) != rM.size1())
                rColumn.resize(rM.size1(), false);
            for (std::size_t i = 0; i < rM.size1(); ++i)
                rColumn[i] = rM(i, J);
        };
        const auto set_column = [](const std::size_t J, DenseMatrixType& rM, const VectorType& rColumn) {
            for (std::size_t i = 0; i < rM.size1(); ++i)
                rM(i, J) = rColumn[i];
        };

        get_column(0, rB, rhs);
        this->InitializeSolutionStep(rA, solution, rhs);

        bool success = true;

        for (std::size_t i_column = 0ul; i_column < n_rhs; ++i_column)
        {
            get_column(i_column, rB, rhs);

            TSparseSpaceType::SetToZero(solution);

            const bool col_success = this->PerformSolutionStep(rA, solution, rhs);
            success = success && col_success;

            set_column(i_column, rX, solution);
        }
    
        this->FinalizeSolutionStep(rA, solution, rhs);
    
        return success;
    
        KRATOS_CATCH("")
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

    static StandardLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenDirectSolver> Factory()
    {
        return StandardLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, EigenDirectSolver>();
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
