//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Denis Demidov
//                   Philipp Bucher
//

// System includes
#include <cstdlib> // for std::getenv

// External includes

// Project includes
#include "input_output/logger.h"
#include "parallel_utilities.h"


namespace Kratos
{

int ParallelUtilities::GetNumThreads()
{
    return GetInstance().msNumThreads;
}

void ParallelUtilities::SetNumThreads(const int NumThreads)
{
#ifdef KRATOS_SMP_NONE
    // do nothing if is shared memory parallelization is disabled
    return;
#endif

    const int num_procs = GetNumProcs();
    KRATOS_WARNING_IF("ParallelUtilities", NumThreads > num_procs) << "The number of requested threads (" << NumThreads << ") exceeds the number of available threads (" << num_procs << ")!" << std::endl;
    GetInstance().msNumThreads = NumThreads;

#if defined(KRATOS_SMP_OPENMP)
    // external libraries included in Kratos still use OpenMP (such as AMGCL)
    // this makes sure that they use the same number of threads as Kratos itself.
    omp_set_num_threads(NumThreads);
#endif
}

int ParallelUtilities::GetNumProcs()
{
#if defined(KRATOS_SMP_OPENMP)
    return omp_get_num_procs();

#elif defined(KRATOS_SMP_CXX11)
    // NOTE: std::thread::hardware_concurrency() can return 0 in some systems!
    unsigned num_procs = std::thread::hardware_concurrency();

    KRATOS_WARNING_IF("ParallelUtilities", num_procs == 0) << "The number of processors cannot be determined correctly on this machine. Please check your setup carefully!" << std::endl;

    return std::max(1, num_procs);

#else
    return 1;
#endif
}

int ParallelUtilities::InitializeNumberOfThreads()
{
#ifdef KRATOS_SMP_NONE
    return 1;

#else
    const char* env_kratos = std::getenv("KRATOS_NUM_THREADS");
    const char* env_omp    = std::getenv("OMP_NUM_THREADS");

    int num_threads;

    if (env_kratos) {
        // "KRATOS_NUM_THREADS" is in the environment
        // giving highest priority to this variable
        num_threads = std::atoi( env_kratos );
        KRATOS_DETAIL("Kernel") << "Using \"KRATOS_NUM_THREADS\" for getting the number of threads: " << num_threads << std::endl;
    } else if (env_omp) {
        // "KRATOS_NUM_THREADS" is not in the environment,
        // checking if "OMP_NUM_THREADS" is
        num_threads = std::atoi( env_omp );
        KRATOS_DETAIL("Kernel") << "Using \"OMP_NUM_THREADS\" for getting the number of threads: " << num_threads << std::endl;
    } else {
        // if neither "KRATOS_NUM_THREADS" not "OMP_NUM_THREADS"
        // is in the environment, then check the C++ thread function
        // NOTE: this can return 0 in some systems!
        num_threads = std::thread::hardware_concurrency();
        KRATOS_DETAIL("Kernel") << "Using \"std::thread::hardware_concurrency\" for getting the number of threads: " << num_threads << std::endl;
    }

    KRATOS_WARNING_IF("ParallelUtilities", num_threads < 1) << "The number of threads could not be determined correctly, which means that Kratos runs in serial / without shared memory parallelism.\nPlease set \"KRATOS_NUM_THREADS\" as environment variable with a value > 0" << std::endl;

    num_threads = std::max(1, num_threads);

#if defined(KRATOS_SMP_OPENMP)
    // external libraries included in Kratos still use OpenMP (such as AMGCL)
    // this makes sure that they use the same number of threads as Kratos itself.
    omp_set_num_threads(num_threads);
#endif

    return num_threads;
#endif
}

ParallelUtilities& ParallelUtilities::GetInstance()
{
    // Using double-checked locking to ensure thread safety in the first creation of the singleton.
    if (!mpInstance) {
        #ifdef KRATOS_SMP_OPENMP
        #pragma omp critical
        if (!mpInstance) {
        #endif
            Create();
            mpInstance->msNumThreads = InitializeNumberOfThreads();
        #ifdef KRATOS_SMP_OPENMP
        }
        #endif
    }

    return *mpInstance;
}

void ParallelUtilities::Create()
{
    static ParallelUtilities parallel_utilities;
    mpInstance = &parallel_utilities;
}

ParallelUtilities* ParallelUtilities::mpInstance = nullptr;
int ParallelUtilities::msNumThreads = 1;

}  // namespace Kratos.
