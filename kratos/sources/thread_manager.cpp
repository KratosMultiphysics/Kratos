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

// System includes
#include <thread>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/thread_manager.h"
#include "input_output/logger.h"

namespace Kratos
{

int ThreadManager::GetNumThreads() const
{
    // KRATOS_ERROR << "Calling the base class ThreadManager::GetNumThreads!" << std::endl;
    KRATOS_WARNING("ThreadManager") << "Calling the base class ThreadManager::GetNumThreads!" << std::endl;
    return 1;
}

/***********************************************************************************/
/***********************************************************************************/

int ThreadManager::GetNumProcs() const
{
    // KRATOS_ERROR << "Calling the base class ThreadManager::GetNumProcs!" << std::endl;
    KRATOS_WARNING("ThreadManager") << "Calling the base class ThreadManager::GetNumProcs!" << std::endl;
    return 1;
}

/***********************************************************************************/
/***********************************************************************************/

void ThreadManager::SetNumThreads(const int NumThreads)
{
    // KRATOS_ERROR << "Calling the base class ThreadManager::SetNumThreads!" << std::endl;
    KRATOS_WARNING("ThreadManager") << "Calling the base class ThreadManager::SetNumThreads!" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

int OMPThreadManager::GetNumThreads() const
{
#ifdef KRATOS_SMP_OPENMP
    int nthreads = omp_get_max_threads();
    KRATOS_DEBUG_ERROR_IF(nthreads <= 0) << "GetNumThreads would devolve nthreads = " << nthreads << " which is not possible" << std::endl;
    return nthreads;
#else
    KRATOS_ERROR << "Calling OMPThreadManager::GetNumThreads when OpenMP is not enabled!" << std::endl;
    return 1;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

int OMPThreadManager::GetNumProcs() const
{
#ifdef KRATOS_SMP_OPENMP
    return omp_get_num_procs();
#else
    KRATOS_ERROR << "Calling OMPThreadManager::GetNumProcs when OpenMP is not enabled!" << std::endl;
    return 1;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

void OMPThreadManager::SetNumThreads(const int NumThreads)
{
#ifdef KRATOS_SMP_OPENMP
    // external libraries included in Kratos still use OpenMP (such as AMGCL)
    // this makes sure that they use the same number of threads as Kratos itself.
    omp_set_num_threads(NumThreads);
#else
    KRATOS_ERROR << "Calling OMPThreadManager::SetNumThreads when OpenMP is not enabled!" << std::endl;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

int CXX11ThreadManager::GetNumThreads() const
{
#if defined(KRATOS_SMP_CXX11)
    int nthreads = std::thread::hardware_concurrency();
    KRATOS_DEBUG_ERROR_IF(nthreads <= 0) << "GetNumThreads would devolve nthreads = " << nthreads << " which is not possible" << std::endl;
    return nthreads;
#else
    KRATOS_ERROR << "Calling CXX11ThreadManager::GetNumThreads when KRATOS_SMP_CXX11 is not enabled!" << std::endl;
    return 1;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

int CXX11ThreadManager::GetNumProcs() const
{
#if defined(KRATOS_SMP_CXX11)
    // NOTE: std::thread::hardware_concurrency() can return 0 in some systems!
    int num_procs = std::thread::hardware_concurrency();

    KRATOS_WARNING_IF("CXX11ThreadManager", num_procs == 0) << "The number of processors cannot be determined correctly on this machine. Please check your setup carefully!" << std::endl;

    return std::max(1, num_procs);
#else
    KRATOS_ERROR << "Calling CXX11ThreadManager::GetNumProcs when KRATOS_SMP_CXX11 is not enabled!" << std::endl;
    return 1;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

void CXX11ThreadManager::SetNumThreads(const int NumThreads)
{
#if defined(KRATOS_SMP_CXX11)
    // NOTE: We do not have a direct way to set the number of threads in C++11
    //       We just store the number of threads in ParallelUtilities and
    //       it is up to the user to use this number when creating threads
    //       (e.g. in IndexPartition or BlockPartition)
#else
    KRATOS_ERROR << "Calling CXX11ThreadManager::SetNumThreads when KRATOS_SMP_CXX11 is not enabled!" << std::endl;
#endif
}

} /// namespace Kratos