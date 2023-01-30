/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Quirin Aumann
*/

#if !defined(KRATOS_EIGEN_DENSE_LLT_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_LLT_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"

namespace Kratos {

template <typename TScalar = double>
class EigenDenseLLTSolver
{
public:
    using Scalar = TScalar;
    using Matrix = Kratos::EigenDynamicMatrix<Scalar>;
    using Vector = Kratos::EigenDynamicVector<Scalar>;

private:
    Eigen::LLT<Matrix> m_solver;

public:
    static std::string Name()
    {
        return "complex_dense_llt";
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

#endif // defined(KRATOS_EIGEN_DENSE_LLT_SOLVER_H_INCLUDED)
