/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Aron Noordam
*/

#pragma once

// External includes
#include <Eigen/Sparse>

// Project includes
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/default_interface.h"
#include "includes/ublas_complex_interface.h"

namespace Kratos {

template <typename TScalar = double>
class EigenSparseLDLTSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = Kratos::EigenSparseMatrix<Scalar>;
    using Vector = Kratos::EigenDynamicVector<Scalar>;

private:
    Eigen::SimplicialLDLT<SparseMatrix> m_solver;

public:
    static std::string Name()
    {
        return "eigen_sparse_ldlt";
    }

    void Initialize(Parameters settings)
    {
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        m_solver.compute(a);

        const bool success = m_solver.info() == Eigen::Success;

        KRATOS_ERROR_IF(!success) << "SimplicialLDLT factorization failed. The matrix might not be symmetric (Hermitian) or might be singular." << std::endl;

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
