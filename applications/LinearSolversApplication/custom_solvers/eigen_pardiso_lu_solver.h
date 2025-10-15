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
#include <Eigen/PardisoSupport>

// Project includes
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "utilities/parallel_utilities.h"

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
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        // Ensure the number of threads in  MKL is the same considered for other operations
        CheckThreadConsistency();

        // Actually compute
        m_solver.compute(a);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    bool Solve(Eigen::Ref<const Vector> b, Eigen::Ref<Vector> x) const
    {
        // Ensure the number of threads in  MKL is the same considered for other operations
        CheckThreadConsistency();
        
        // Actually solve
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

private:

    /**
     * @brief Checks and adjusts the number of threads used by the MKL library for consistency.
     * @details This method compares the maximum number of threads configured in MKL (Intel Math Kernel Library) with the number of threads currently being used by the application (obtained via `ParallelUtilities::GetNumThreads()`).
     * If MKL's thread count is greater than the application's thread count, a warning is issued, and MKL's thread count is reduced to match the application's thread count via `mkl_set_num_threads()`.
     * This ensures that MKL operations (potentially used by the solver) do not utilize more threads than allowed by the overall parallel configuration, preventing potential oversubscription or inconsistency.
     */
    static void CheckThreadConsistency()
    {
        const int number_of_threads_mkl = mkl_get_max_threads();
        const int number_of_threads_used = ParallelUtilities::GetNumThreads();
        if (number_of_threads_mkl > number_of_threads_used) {
            KRATOS_WARNING("EigenPardisoLUSolver") << "Setting the number of threads in MKL to adapt to ParallelUtilities::GetNumThreads(): " << number_of_threads_used << " instead of " << number_of_threads_mkl << std::endl;
            mkl_set_num_threads(number_of_threads_used);
        }
    }
};

} // namespace Kratos
