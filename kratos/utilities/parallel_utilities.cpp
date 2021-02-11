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
#include <algorithm>
#include <cstdlib> // for std::getenv

// External includes

// Project includes
#include "parallel_utilities.h"
#include "input_output/logger.h"
#include "includes/lock_object.h"


namespace Kratos
{

int ParallelUtilities::GetNumThreads()
{
    return GetNumberOfThreads();
}

void ParallelUtilities::SetNumThreads(const int NumThreads)
{
#ifdef KRATOS_SMP_NONE
    // do nothing if is shared memory parallelization is disabled
    return;
#endif

    const int num_procs = GetNumProcs();
    KRATOS_WARNING_IF("ParallelUtilities", NumThreads > num_procs) << "The number of requested threads (" << NumThreads << ") exceeds the number of available threads (" << num_procs << ")!" << std::endl;
    GetNumberOfThreads() = NumThreads;

#ifdef KRATOS_SMP_OPENMP
    // external libraries included in Kratos still use OpenMP (such as AMGCL)
    // this makes sure that they use the same number of threads as Kratos itself.
    omp_set_num_threads(NumThreads);
#endif
}

int ParallelUtilities::GetNumProcs()
{
#ifdef KRATOS_SMP_OPENMP
    return omp_get_num_procs();

#elif defined(KRATOS_SMP_CXX11)
    // NOTE: std::thread::hardware_concurrency() can return 0 in some systems!
    int num_procs = std::thread::hardware_concurrency();

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
    } else if (env_omp) {
        // "KRATOS_NUM_THREADS" is not in the environment,
        // checking if "OMP_NUM_THREADS" is
        num_threads = std::atoi( env_omp );
    } else {
        // if neither "KRATOS_NUM_THREADS" not "OMP_NUM_THREADS"
        // is in the environment, then check the C++ thread function
        // NOTE: this can return 0 in some systems!
        num_threads = std::thread::hardware_concurrency();
    }

    num_threads = std::max(1, num_threads);

#ifdef KRATOS_SMP_OPENMP
    // external libraries included in Kratos still use OpenMP (such as AMGCL)
    // this makes sure that they use the same number of threads as Kratos itself.
    omp_set_num_threads(num_threads);
#endif

    return num_threads;
#endif
}

int& ParallelUtilities::GetNumberOfThreads()
{
    if (!mspNumThreads) {
        LockObject lock;
        lock.SetLock();
        if (!mspNumThreads) {
            static int num_threads;
            num_threads = InitializeNumberOfThreads();
            mspNumThreads = &num_threads;
        }
        lock.UnSetLock();
    }

    return *mspNumThreads;
}

int* ParallelUtilities::mspNumThreads = nullptr;

}  // namespace Kratos.
