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

namespace Kratos
{

/**
 * @class MKLUtilities
 * @brief Utility class for managing Intel MKL (Math Kernel Library) threading configuration
 * @details This class provides static methods and enumerations to handle thread management for Intel MKL operations
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(LINEARSOLVERS_APPLICATION) MKLUtilities
{
public:
    /**
     * @brief Enumeration for MKL thread setting options in Pardiso solver
     * @details This enum defines different threading strategies that can be used with the Intel MKL Pardiso solver to control parallel execution behavior.
     */
    enum class MKLThreadSetting
    {
        Minimal = -2,    /// Use minimal number of threads for computation
        Consistent = -1, /// Use consistent thread count across solver operations
        Do_nothing = 0,  /// Do not change the number of threads
        Manual = 1       /// Manually specify the number of threads to use
    };

    static int GetNumMKLThreads();

    static void SetNumMKLThreads(const int NumThreads);

    /**
     * @brief Checks if the MKL thread number has been previously validated or set.
     * @details This method determines whether the MKL thread configuration has already been checked or configured during the current execution. It serves as a guard to prevent redundant thread setting operations.
     * @param NumberOfMKLThreads The MKL thread setting value to check against
     * @return true if the thread number has been previously checked/set, false if it needs to be configured
     * @note This is typically used internally by SetMKLThreadCount() to ensure thread configuration happens only once per execution
     */
    static bool CheckThreadNumber(const int NumberOfMKLThreads);

    /**
     * @brief Sets the number of threads for MKL (Intel Math Kernel Library) operations.
     * @details This method configures the MKL thread count based on the specified thread setting policy.
     * Thread count is only set once per execution (controlled by CheckThreadNumber()).
     * The thread count is determined by the following policies:
     * - Positive value (Manual): Uses the explicitly specified number of MKL threads
     * - MKLThreadSetting::Minimal (-2): Uses the minimum between available MKL threads and currently used threads
     * - MKLThreadSetting::Consistent (-1): Uses the number of threads from ParallelUtilities
     * - MKLThreadSetting::Do_nothing (0): No thread count modification is performed
     * - Invalid setting: Throws a KRATOS_ERROR
     * @param NumberOfMKLThreads The thread setting value (positive integer or MKLThreadSetting enum value)
     * @note This method modifies the global MKL thread count using mkl_set_num_threads()
     * @note Thread count configuration is performed only once per execution to avoid conflicts
     */
    static void SetMKLThreadCount(const int NumberOfMKLThreads);

};

} // namespace Kratos
