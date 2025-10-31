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

// Project includes
#include "custom_utilities/mkl_utilities.h"

namespace Kratos
{

bool MKLUtilities::CheckThreadNumber(const int NumberOfMKLThreads)
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

/***********************************************************************************/
/***********************************************************************************/

void MKLUtilities::SetMKLThreadCount(const int NumberOfMKLThreads)
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

} // namespace Kratos
