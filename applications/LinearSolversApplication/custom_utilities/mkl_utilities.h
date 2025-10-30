//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// External includes
#include <mkl.h>

namespace Kratos
{

class MKLUtilities
{
public:
    /**
     * @brief Enumeration for MKL thread setting options in Pardiso solver
     * @details This enum defines different threading strategies that can be used with the Intel MKL Pardiso solver to control parallel execution behavior.
     */
    enum class MKLThreadSetting
    {
        Minimal = -1,   /// Use minimal number of threads for computation
        Consistent = 0, /// Use consistent thread count across solver operations
        Manual = 1      /// Manually specify the number of threads to use
    };

    /**
     * @brief Checks the consistency between MKL thread configuration and the application's thread count.
     * @details This method verifies that the MKL Pardiso solver's thread configuration aligns with the application's threading policy based on the configured MKL thread setting:
     * Manual mode (NumberOfMKLThreads > 0):
     * - Checks if MKL's current thread count matches the manually specified number
     * - Issues warning if there's a mismatch
     * Minimal mode (MKLThreadSetting::Minimal):
     * - Ensures MKL uses the minimum between its current setting and application threads
     * - Issues warning if MKL thread count exceeds this minimum
     * Consistent mode (MKLThreadSetting::Consistent):
     * - Ensures MKL thread count doesn't exceed the application's thread count
     * - Issues warning if MKL is configured to use more threads than available to the application
     * @note This method only performs validation and issues warnings; it does NOT modify the MKL thread count. Use EnsureMKLThreadConsistency() or SetMKLThreadCount() to apply corrections.
     * @return true if the thread configuration is consistent, false if adjustment is needed
     */
    static bool CheckThreadNumber(const int NumberOfMKLThreads)
    {
        const int number_of_threads_mkl = mkl_get_max_threads();
        if (NumberOfMKLThreads > 0) { // Manual setting
            if (number_of_threads_mkl != NumberOfMKLThreads) {
                KRATOS_WARNING("EigenPardisoLUSolver") << "The number of threads in MKL is: " << NumberOfMKLThreads << " instead of " << number_of_threads_mkl << std::endl;
                return false;
            }
            return true;
        } else { // NumberOfMKLThreads <= 0, for Minimal or Consistent
            const int number_of_threads_used = ParallelUtilities::GetNumThreads();
            if (static_cast<int>(MKLThreadSetting::Minimal) == NumberOfMKLThreads) {
                const int min_threads = std::min(number_of_threads_mkl, number_of_threads_used);
                if (number_of_threads_mkl != min_threads) {
                    KRATOS_WARNING("EigenPardisoLUSolver") << "The number of threads in MKL is: " << number_of_threads_mkl << " instead of minimal: " << min_threads << std::endl;
                    return false;
                }
                return true;
            } else if (static_cast<int>(MKLThreadSetting::Consistent) == NumberOfMKLThreads) {
                if (number_of_threads_mkl > number_of_threads_used) {
                    KRATOS_WARNING("EigenPardisoLUSolver") << "The number of threads in MKL is: " << number_of_threads_mkl << " which is greater than the number of threads used by the application: " << number_of_threads_used << std::endl;
                    return false;
                }
                return true;
            } else {
                KRATOS_ERROR << "Invalid MKL thread setting: " << NumberOfMKLThreads << std::endl;
            }
        }
        return true;
    }

    /**
     * @brief Sets the number of threads for MKL (Intel Math Kernel Library) operations.
     * @details This method configures the MKL thread count based on the current thread setting policy.
     * It only executes if the thread number hasn't been previously checked/set.
     * The thread count is determined by the following policies:
     * - Positive value: Uses the explicitly set number of MKL threads
     * - MKLThreadSetting::Minimal: Uses the minimum between available MKL threads and currently used threads
     * - MKLThreadSetting::Consistent: Uses the number of threads from ParallelUtilities
     * - Invalid setting: Throws a KRATOS_ERROR
     * @note This method modifies the global MKL thread count using mkl_set_num_threads()
     * @note Thread count is only set once per execution (controlled by CheckThreadNumber())
     */
    static void SetMKLThreadCount(const int NumberOfMKLThreads)
    {
        if (!CheckThreadNumber(NumberOfMKLThreads)) {
            int number_of_threads_used;
            if (NumberOfMKLThreads > 0) {
                number_of_threads_used = NumberOfMKLThreads;
            } else if (static_cast<int>(MKLThreadSetting::Minimal) == NumberOfMKLThreads) {
                const int number_of_threads_mkl = mkl_get_max_threads();
                number_of_threads_used = std::min(number_of_threads_mkl, number_of_threads_used);
            } else if (static_cast<int>(MKLThreadSetting::Consistent) == NumberOfMKLThreads) {
                number_of_threads_used = ParallelUtilities::GetNumThreads();
            } else {
                KRATOS_ERROR << "Invalid MKL thread setting: " << NumberOfMKLThreads << std::endl;
            }
            mkl_set_num_threads(number_of_threads_used);
        }
    }
};

} // namespace Kratos
