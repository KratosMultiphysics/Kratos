/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _ | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#pragma once

// System includes
#include <optional>

// External includes
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

// Project includes
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "custom_utilities/mkl_utilities.h"

namespace Kratos {

template <typename TScalar = double>
class EigenPardisoLUSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = Kratos::EigenSparseMatrix<Scalar>;
    using Vector = Kratos::EigenDynamicVector<Scalar>;

private:
    Eigen::PardisoLU<SparseMatrix> m_solver;

public:
    static std::string Name()
    {
        return "eigen_pardiso_lu";
    }

    void Initialize(Parameters settings)
    {
        // Compute the number of MKL threads
        mNumberOfMKLThreads = MKLUtilities::ComputeMKLThreadCount(settings);
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        const int previous_threads = MKLUtilities::GetNumThreads();
        if (mNumberOfMKLThreads) MKLUtilities::SetNumThreads(*mNumberOfMKLThreads);
        m_solver.compute(a);
        if (mNumberOfMKLThreads) MKLUtilities::SetNumThreads(previous_threads);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    bool Solve(Eigen::Ref<const Vector> b, Eigen::Ref<Vector> x) const
    {
        const int previous_threads = MKLUtilities::GetNumThreads();
        if (mNumberOfMKLThreads) MKLUtilities::SetNumThreads(*mNumberOfMKLThreads);
        x = m_solver.solve(b);
        if (mNumberOfMKLThreads) MKLUtilities::SetNumThreads(previous_threads);

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

private:
    ///@name Private Member Variables
    ///@{

    std::optional<int> mNumberOfMKLThreads = std::nullopt; /// The number of MKL threads to be used

    ///@}
};

} // namespace Kratos
