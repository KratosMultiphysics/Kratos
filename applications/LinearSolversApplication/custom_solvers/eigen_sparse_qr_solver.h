/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_EIGEN_SPARSE_QR_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_SPARSE_QR_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

namespace Kratos {

template <typename TScalar = double>
class EigenSparseQRSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = Kratos::EigenSystemSparseMatrix<Scalar>;
    using Vector = Kratos::EigenDynamicVector<Scalar>;

private:
    // The fill-reducing ordering must use the same index type as the matrix
    // (int for the uBLAS backend, the system index for the Eigen backend).
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<typename SparseMatrix::StorageIndex>> m_solver;

public:
    static std::string Name()
    {
        return "eigen_sparse_qr";
    }

    void Initialize(Parameters settings)
    {
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        m_solver.compute(a);

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
        return m_solver.lastErrorMessage();
    }
};

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_SPARSE_QR_SOLVER_H_INCLUDED)
