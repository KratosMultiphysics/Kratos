/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _ | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
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
        // Ensure the number of threads in  MKL is the same considered for other operations
        EnsureMKLThreadConsistency();
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        // Check thread consistency
        CheckThreadConsistency();

        // Actually compute
        m_solver.compute(a);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    bool Solve(Eigen::Ref<const Vector> b, Eigen::Ref<Vector> x) const
    {
        // Check thread consistency
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
     * @brief Checks the consistency between the maximum number of threads configured in MKL and the application's thread count.
     * @details This method compares the maximum number of threads configured in MKL (mkl_get_max_threads()) with the number of threads currently used by the application (ParallelUtilities::GetNumThreads()).
     * - If MKL threads > Application threads: A warning is issued, indicating that MKL is configured to use more threads than the application allows. The method returns false, signaling to the caller that the MKL thread count needs to be adjusted (reduced).
     * - If MKL threads <= Application threads: The thread count is considered consistent. No action is taken, and the method returns true.
     * @note This method only performs the check and warning; it does NOT modify the MKL thread count.** The caller is responsible for applying the necessary adjustment if the check returns false.
     * @return true if MKL's thread count is consistent (less than or equal to the application's) or false if it is inconsistent (MKL thread count is greater).
     */
    static bool CheckThreadConsistency()
    {
        const int number_of_threads_mkl = mkl_get_max_threads();
        const int number_of_threads_used = ParallelUtilities::GetNumThreads();
        if (number_of_threads_mkl > number_of_threads_used) {
            KRATOS_WARNING("EigenPardisoLUSolver") << "Setting the number of threads in MKL to adapt to ParallelUtilities::GetNumThreads(): " << number_of_threads_used << " instead of " << number_of_threads_mkl << std::endl;
            return false;
        }
        return true;
    }
    
    /**
    * @brief Ensures that MKL's thread count does not exceed the application's configured thread count.
    * @details Calls CheckThreadConsistency(). If the check returns false (indicating MKL threads > Application threads), it reduces MKL's thread count to match the application's count via mkl_set_num_threads().
    * This method effectively performs the correction identified by CheckThreadConsistency().
    */
    static void EnsureMKLThreadConsistency()
    {
        if (!CheckThreadConsistency()) {
            const int number_of_threads_used = ParallelUtilities::GetNumThreads();
            mkl_set_num_threads(number_of_threads_used);
        }
    }
};

} // namespace Kratos
