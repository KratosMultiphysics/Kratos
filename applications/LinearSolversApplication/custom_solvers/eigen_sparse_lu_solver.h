/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#pragma once

// External includes
#include <Eigen/Sparse>

// Project includes
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

namespace Kratos {

template <typename TScalar = double>
class EigenSparseLUSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = Kratos::EigenSystemSparseMatrix<Scalar>;
    using Vector = Kratos::EigenDynamicVector<Scalar>;

private:
    // NOTE: Eigen::SparseLU requires its matrix type to be column-major (see Eigen's
    // SparseLU::factorize(), which assumes matrix.outerIndexPtr() are column pointers
    // when applying the column permutation). Kratos' EigenSparseMatrix is row-major
    // (to match Kratos'/uBLAS' row-major convention), so we convert to a column-major
    // copy before handing it to the solver, rather than instantiating SparseLU with a
    // row-major matrix type (which silently corrupts its internal state and crashes).
    using ColMajorSparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, int>;

    Eigen::SparseLU<ColMajorSparseMatrix> m_solver;

public:
    static std::string Name()
    {
        return "eigen_sparse_lu";
    }

    void Initialize(Parameters settings)
    {
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        const ColMajorSparseMatrix a_col = a;

        m_solver.compute(a_col);

        const bool success = m_solver.info() == Eigen::Success;

        KRATOS_ERROR_IF(!success) << m_solver.lastErrorMessage() << std::endl;

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
        return m_solver.lastErrorMessage();
    }
};

} // namespace Kratos