/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Quirin Aumann
*/

#if !defined(KRATOS_EIGEN_DENSE_COL_PIV_HOUSEHOLDER_QR_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_COL_PIV_HOUSEHOLDER_QR_SOLVER_H_INCLUDED

// External includes
#include <Eigen/QR>

// Project includes
#include "includes/define.h"

namespace Kratos {

template <typename TScalar = double>
class EigenDenseColPivHouseholderQRSolver
{
public:
    using Scalar = TScalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

private:
    Eigen::ColPivHouseholderQR<Matrix> m_solver;

public:
    static std::string Name()
    {
        return "complex_dense_col_piv_householder_qr";
    }

    void Initialize(Parameters settings)
    {
    }

    bool Compute(Eigen::Map<Matrix> a)
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

    bool SolveMultiple(Eigen::Ref<const Matrix> b, Eigen::Ref<Matrix> x) const
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

#endif // defined(KRATOS_EIGEN_DENSE_COL_PIV_HOUSEHOLDER_QR_SOLVER_H_INCLUDED)
