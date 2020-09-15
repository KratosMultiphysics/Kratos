//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


// System includes
#include <cstdlib>
#include <thread>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "input_output/logger.h"
#include "utilities/parallel_helpers.h"


namespace Kratos {

int ParallelHelpers::GetNumThreads()
{
#if defined(KRATOS_SMP_OPENMP)
    return omp_get_max_threads();

#elif defined(KRATOS_SMP_CXX11)
    const char* env_kratos = std::getenv("KRATOS_NUM_THREADS");
    const char* env_omp    = std::getenv("OMP_NUM_THREADS");

    int num_threads = 1;

    if (env_kratos) {
        // "KRATOS_NUM_THREADS" is in the environment
        // giving highest priority to this variable
        num_threads = std::atoi( env_kratos );
        KRATOS_DETAIL("ParallelHelpers") << "Using \"KRATOS_NUM_THREADS\" for \"GetNumThreads\": " << num_threads << std::endl;
    } else if (env_omp) {
        // "KRATOS_NUM_THREADS" is not in the environment,
        // checking if "OMP_NUM_THREADS" is
        num_threads = std::atoi( env_omp );
        KRATOS_DETAIL("ParallelHelpers") << "Using \"OMP_NUM_THREADS\" for \"GetNumThreads\": " << num_threads << std::endl;
    } else {
        // if neither "KRATOS_NUM_THREADS" not "OMP_NUM_THREADS"
        // is in the environment, then check the C++ thread function
        // NOTE: this can return 0 in some systems!
        num_threads = std::thread::hardware_concurrency();
        KRATOS_DETAIL("ParallelHelpers") << "Using \"std::thread::hardware_concurrency\" for \"GetNumThreads\": " << num_threads << std::endl;
    }

    return std::max(1, num_threads);

#else
    return 1;

#endif
}


}  // namespace Kratos.


