//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Denis Demidov
//                   Philipp Bucher
//

// System includes
#include <algorithm>
#include <cstdlib> // for std::getenv
#include <sstream> // for std::stringstream

// External includes

// Project includes
#include "parallel_utilities.h"
#include "input_output/logger.h"
#include "includes/lock_object.h"

namespace Kratos {

namespace {
    std::once_flag flag_once;
}

#ifdef KRATOS_SMP_OPENMP
std::vector<ThreadManager::Pointer> ParallelUtilities::msThreadManagers = { Kratos::make_shared<OMPThreadManager>() };
#elif defined(KRATOS_SMP_CXX11)
std::vector<ThreadManager::Pointer> ParallelUtilities::msThreadManagers = { Kratos::make_shared<CXX11ThreadManager>() };
#else // KRATOS_SMP_NONE
std::vector<ThreadManager::Pointer> ParallelUtilities::msThreadManagers = { Kratos::make_shared<ThreadManager>() };
#endif

int ParallelUtilities::GetNumThreads()
{
    // Get the number of threads available in each ThreadManager
    std::vector<int> num_threads;
    num_threads.reserve(msThreadManagers.size());
    for (const auto& r_manager : msThreadManagers) {
        num_threads.push_back(r_manager->GetNumThreads());
    }

    // Throw a warning if the number of threads is not the same in all ThreadManagers
    if (!std::all_of(num_threads.begin(), num_threads.end(), [&num_threads](const int value){ return value == num_threads[0]; })) {
        std::stringstream err_msg;
        err_msg << "The number of threads is not the same in all ThreadManagers:\n";
        for (std::size_t i=0; i<num_threads.size(); ++i) {
            err_msg << "ThreadManager " << msThreadManagers[i]->Info() << ": " << num_threads[i] << " threads.\n";
        }
        KRATOS_WARNING("ParallelUtilities") << err_msg.str();
    }

    // Return the minimum number of threads available
    return *std::min_element(num_threads.begin(), num_threads.end());
}

void ParallelUtilities::SetNumThreads(const int NumThreads)
{
    KRATOS_ERROR_IF(NumThreads <= 0) << "Attempting to set NumThreads to <= 0. This is not allowed" << std::endl;

#ifdef KRATOS_SMP_NONE
    // do nothing if is shared memory parallelization is disabled
    return;
#endif

    const int num_procs = GetNumProcs();
    KRATOS_WARNING_IF("ParallelUtilities", NumThreads > num_procs) << "The number of requested threads (" << NumThreads << ") exceeds the number of available threads (" << num_procs << ")!" << std::endl;
    GetNumberOfThreads() = NumThreads;

    // Set the number of threads in the ThreadManagers
    for (const auto& r_manager : msThreadManagers) {
        r_manager->SetNumThreads(NumThreads);
    }
}

int ParallelUtilities::GetNumProcs()
{
    // Get the number of processors available in each ThreadManager
    std::vector<int> num_procs;
    num_procs.reserve(msThreadManagers.size());
    for (const auto& r_manager : msThreadManagers) {
        num_procs.push_back(r_manager->GetNumProcs());
    }

    // Throw a warning if the number of processors is not the same in all ThreadManagers
    if (!std::all_of(num_procs.begin(), num_procs.end(), [&num_procs](const int value){ return value == num_procs[0]; })) {
        std::stringstream err_msg;
        err_msg << "The number of processors is not the same in all ThreadManagers:\n";
        for (std::size_t i=0; i<num_procs.size(); ++i) {
            err_msg << "ThreadManager " << msThreadManagers[i]->Info() << ": " << num_procs[i] << " processors.\n";
        }
        KRATOS_WARNING("ParallelUtilities") << err_msg.str();
    }

    // Return the minimum number of processors available
    return *std::min_element(num_procs.begin(), num_procs.end());
}

int ParallelUtilities::InitializeNumberOfThreads()
{
#ifdef KRATOS_SMP_NONE
    return 1;
#else
    const char* env_kratos = std::getenv("KRATOS_NUM_THREADS");
    const char* env_omp    = std::getenv("OMP_NUM_THREADS");

    int num_threads;

    // Determine the number of threads to use
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

    // Ensure at least one thread is available
    num_threads = std::max(1, num_threads);

    // Initialize the number of threads in the ThreadManagers
    for (const auto& r_manager : msThreadManagers) {
        r_manager->SetNumThreads(num_threads);
    }

    return num_threads;
#endif
}

LockObject& ParallelUtilities::GetGlobalLock()
{
    if (!mspGlobalLock) {
        std::call_once(flag_once, [](){
            static LockObject global_lock;
            mspGlobalLock = &global_lock;
        });
    }

    return *mspGlobalLock;
}

void ParallelUtilities::AddThreadManager(ThreadManager::Pointer pThreadManager)
{
    msThreadManagers.push_back(pThreadManager);
}

int& ParallelUtilities::GetNumberOfThreads()
{
    if (!mspNumThreads) {
        KRATOS_CRITICAL_SECTION
        if (!mspNumThreads) {
            static int num_threads;
            num_threads = InitializeNumberOfThreads();
            mspNumThreads = &num_threads;
        }
    }

    return *mspNumThreads;
}

LockObject* ParallelUtilities::mspGlobalLock = nullptr;
int* ParallelUtilities::mspNumThreads = nullptr;

}  // namespace Kratos.
