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

#if !defined(KRATOS_EIGEN_SPARSE_CG_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_SPARSE_CG_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

namespace Kratos {

template <typename TScalar = double>
class EigenSparseCGSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor, int>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

private:
    SparseMatrix m_a;
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> m_solver;

public:
    static std::string Name()
    {
        return "eigen_sparse_cg";
    }

    void Initialize(Parameters settings)
    {
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        m_a = a;

        m_solver.compute(m_a);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    bool Solve(Eigen::Ref<const Vector> b, Eigen::Ref<Vector> x) const
    {
        x = m_solver.solve(b);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "EigenDirectSolver <" << Name() << "> finished.";
    }

    std::string GetSolverErrorMessages() const
    {
        return "No additional information";
    }
};

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_SPARSE_CG_SOLVER_H_INCLUDED)
