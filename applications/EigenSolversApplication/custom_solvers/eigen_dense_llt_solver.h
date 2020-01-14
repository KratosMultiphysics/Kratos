/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Quirin Aumann
*/

#if !defined(KRATOS_EIGEN_DENSE_LLT_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_LLT_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos {

template <typename TScalar = double>
class EigenDenseLLTSolver
{
public:
    using Scalar = TScalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

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
