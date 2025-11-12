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

// System includes
#include <optional>

// Project includes
#include "includes/kratos_parameters.h"

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

    /**
     * @brief Retrieves the current number of threads configured for Intel MKL operations.
     * @details This static method queries the Intel Math Kernel Library (MKL) to obtain the number of threads currently set for parallel computations.
     * This value reflects the thread configuration that will be used for subsequent MKL operations in the current process.
     * @return The number of threads currently set for MKL computations.
     */
    static int GetNumThreads();

    /**
     * @brief Sets the number of threads to be used by Intel MKL library operations.
     * @details This static method configures the maximum number of threads that Intel MKL (Math Kernel Library) will use for parallel computations. This setting affects all subsequent MKL operations in the current process.
     * @param NumThreads The number of threads to be used by MKL. Typically set to the number of available CPU cores or a fraction thereof depending on the desired parallelization strategy.
     */
    static void SetNumThreads(const unsigned int NumThreads);

    /**
     * @brief Checks if the MKL thread number has been previously validated or set.
     * @details This method determines whether the MKL thread configuration has already been checked or configured during the current execution. It serves as a guard to prevent redundant thread setting operations.
     * @param NumberOfMKLThreads The MKL thread setting value to check against
     * @return true if the thread number has been previously checked/set, false if it needs to be configured
     * @note This is typically used internally by ComputeMKLThreadCount() to ensure thread configuration happens only once per execution
     */
    static bool CheckThreadNumber(const int NumberOfMKLThreads);

    /**
     * @brief Computes and configures the optimal number of threads for MKL operations. Version with Parameters input.
     * @details This method determines and sets the appropriate MKL thread count based on the specified threading policy.
     * The configuration is applied only once per execution (guarded by CheckThreadNumber()) to prevent conflicts.
     * Threading policies:
     * - Positive integer (Manual mode): Explicitly sets the specified number of MKL threads
     * - MKLThreadSetting::Minimal (-2): Sets threads to minimum of available MKL threads and current thread count
     * - MKLThreadSetting::Consistent (-1): Aligns MKL threads with ParallelUtilities thread count
     * - MKLThreadSetting::Do_nothing (0): Preserves current MKL thread configuration
     * - Invalid values: Triggers KRATOS_ERROR exception
     * @param Settings The settings including the options.
     * @return The final number of threads configured for MKL operations.
     */
    static std::optional<int> ComputeMKLThreadCount(Parameters Settings);

    /**
     * @brief Computes and configures the optimal number of threads for MKL operations
     * @details This method determines and sets the appropriate MKL thread count based on the specified threading policy.
     * The configuration is applied only once per execution (guarded by CheckThreadNumber()) to prevent conflicts.
     * Threading policies:
     * - Positive integer (Manual mode): Explicitly sets the specified number of MKL threads
     * - MKLThreadSetting::Minimal (-2): Sets threads to minimum of available MKL threads and current thread count
     * - MKLThreadSetting::Consistent (-1): Aligns MKL threads with ParallelUtilities thread count
     * - MKLThreadSetting::Do_nothing (0): Preserves current MKL thread configuration
     * - Invalid values: Triggers KRATOS_ERROR exception
     * @param NumberOfMKLThreads Thread configuration parameter (positive integer or MKLThreadSetting enum cast to int).
     * @return The final number of threads configured for MKL operations.
     */
    static std::optional<int> ComputeMKLThreadCount(const int NumberOfMKLThreads);

};

} // namespace Kratos
