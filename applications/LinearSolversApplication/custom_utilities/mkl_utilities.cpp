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

// External includes
#include <mkl.h>

// Project includes
#include "input_output/logger.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/mkl_utilities.h"

namespace Kratos
{

int MKLUtilities::GetNumThreads()
{
    return mkl_get_max_threads();
}

/***********************************************************************************/
/***********************************************************************************/

void MKLUtilities::SetNumThreads(const unsigned int NumThreads)
{
    mkl_set_num_threads(NumThreads);
}

/***********************************************************************************/
/***********************************************************************************/

bool MKLUtilities::CheckThreadNumber(const int NumberOfMKLThreads)
{
    // Do nothing case
    if (NumberOfMKLThreads == static_cast<int>(MKLThreadSetting::Do_nothing)) {
        return true;
    }

    const int number_of_threads_mkl = GetNumThreads();
    if (NumberOfMKLThreads > 0) { // Manual setting
        if (number_of_threads_mkl != NumberOfMKLThreads) {
            KRATOS_WARNING("MKLUtilities") << "The number of threads in MKL is: " << NumberOfMKLThreads << " instead of " << number_of_threads_mkl << std::endl;
            return false;
        }
        return true;
    } else { // NumberOfMKLThreads <= 0, for Minimal or Consistent
        const int number_of_threads_used = ParallelUtilities::GetNumThreads();
        if (static_cast<int>(MKLThreadSetting::Minimal) == NumberOfMKLThreads) {
            const int min_threads = std::min(number_of_threads_mkl, number_of_threads_used);
            if (number_of_threads_mkl != min_threads) {
                KRATOS_WARNING("MKLUtilities") << "The number of threads in MKL is: " << number_of_threads_mkl << " instead of minimal: " << min_threads << std::endl;
                return false;
            }
            return true;
        } else if (static_cast<int>(MKLThreadSetting::Consistent) == NumberOfMKLThreads) {
            if (number_of_threads_mkl > number_of_threads_used) {
                KRATOS_WARNING("MKLUtilities") << "The number of threads in MKL is: " << number_of_threads_mkl << " which is greater than the number of threads used by the application: " << number_of_threads_used << std::endl;
                return false;
            }
            return true;
        } else {
            KRATOS_ERROR << "Invalid MKL thread setting: " << NumberOfMKLThreads << std::endl;
        }
    }
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

std::optional<int> MKLUtilities::ComputeMKLThreadCount(Parameters Settings)
{
    // Set default parameters
    if (!Settings.Has("num_threads_mkl")) {
        Settings.AddEmptyValue("num_threads_mkl").SetString("do_nothing"); // Default to 0 (do nothing)
    }

    // Configure number of threads for MKL Pardiso solver
    int number_of_mkl_threads = 0;
    if (Settings["num_threads_mkl"].IsNumber()) {
        number_of_mkl_threads = Settings["num_threads_mkl"].GetInt();
    } else if (Settings["num_threads_mkl"].GetString() == "minimal") {
        number_of_mkl_threads = static_cast<int>(MKLUtilities::MKLThreadSetting::Minimal);
    } else if (Settings["num_threads_mkl"].GetString() == "consistent") {
        number_of_mkl_threads = static_cast<int>(MKLUtilities::MKLThreadSetting::Consistent);
    } else if (Settings["num_threads_mkl"].GetString() == "do_nothing") {
        number_of_mkl_threads = static_cast<int>(MKLUtilities::MKLThreadSetting::Do_nothing);
    } else {
        KRATOS_ERROR << "Invalid value for 'num_threads_mkl': " << Settings["num_threads_mkl"].GetString() << ". Accepted values are 'minimal', 'consistent', or an integer." << std::endl;
    }

    // Ensure the number of threads in MKL is the same considered for other operations
    return ComputeMKLThreadCount(number_of_mkl_threads);
}

/***********************************************************************************/
/***********************************************************************************/

std::optional<int> MKLUtilities::ComputeMKLThreadCount(const int NumberOfMKLThreads)
{
    // Check first if it is needed to set the number of threads
    if (!CheckThreadNumber(NumberOfMKLThreads)) {
        int number_of_threads_used;
        if (NumberOfMKLThreads > 0) {
            number_of_threads_used = NumberOfMKLThreads;
        } else if (static_cast<int>(MKLThreadSetting::Minimal) == NumberOfMKLThreads) {
            const int number_of_threads_mkl = GetNumThreads();
            number_of_threads_used = std::min(number_of_threads_mkl, number_of_threads_used);
        } else if (static_cast<int>(MKLThreadSetting::Consistent) == NumberOfMKLThreads) {
            number_of_threads_used = ParallelUtilities::GetNumThreads();
        } else {
            KRATOS_ERROR << "Invalid MKL thread setting: " << NumberOfMKLThreads << std::endl;
        }
        return number_of_threads_used;
    }
    return std::nullopt;
}

} // namespace Kratos
