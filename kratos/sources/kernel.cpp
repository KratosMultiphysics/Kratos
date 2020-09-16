//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes
#include <iostream>
#include <cstdlib>
#include <thread>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/kernel.h"
#include "includes/kratos_version.h"
#include "includes/data_communicator.h"
#include "includes/parallel_environment.h"
#include "input_output/logger.h"
#include "utilities/openmp_utils.h"

namespace Kratos {

namespace {

int DetermineNumberOfThreads()
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
        KRATOS_DETAIL("Kernel") << "Using \"KRATOS_NUM_THREADS\" for \"GetNumThreads\": " << num_threads << std::endl;
    } else if (env_omp) {
        // "KRATOS_NUM_THREADS" is not in the environment,
        // checking if "OMP_NUM_THREADS" is
        num_threads = std::atoi( env_omp );
        KRATOS_DETAIL("Kernel") << "Using \"OMP_NUM_THREADS\" for \"GetNumThreads\": " << num_threads << std::endl;
    } else {
        // if neither "KRATOS_NUM_THREADS" not "OMP_NUM_THREADS"
        // is in the environment, then check the C++ thread function
        // NOTE: this can return 0 in some systems!
        num_threads = std::thread::hardware_concurrency();
        KRATOS_DETAIL("Kernel") << "Using \"std::thread::hardware_concurrency\" for \"GetNumThreads\": " << num_threads << std::endl;
    }

    return std::max(1, num_threads);

#else
    return 1;

#endif
}

}


Kernel::Kernel(bool IsDistributedRun) : mpKratosCoreApplication(Kratos::make_shared<KratosApplication>(
                std::string("KratosMultiphysics"))) {
    mIsDistributedRun = IsDistributedRun;
    KRATOS_INFO("") << " |  /           |\n"
                    << " ' /   __| _` | __|  _ \\   __|\n"
                    << " . \\  |   (   | |   (   |\\__ \\\n"
                    << "_|\\_\\_|  \\__,_|\\__|\\___/ ____/\n"
                    << "           Multi-Physics " << GetVersionString() << std::endl;

    PrintParallelismSupportInfo();

    if (!IsImported("KratosMultiphysics")) {
        this->ImportApplication(mpKratosCoreApplication);
    }

    mNumThreads = DetermineNumberOfThreads();
}

std::unordered_set<std::string> &Kernel::GetApplicationsList() {
  static std::unordered_set<std::string> application_list;
  return application_list;
}

bool Kernel::IsImported(std::string ApplicationName) const {
    if (GetApplicationsList().find(ApplicationName) !=
        GetApplicationsList().end())
        return true;
    else
        return false;
}

bool Kernel::IsDistributedRun() {
    return mIsDistributedRun;
}

int Kernel::GetNumThreads()
{
    return mNumThreads;
}

void Kernel::SetNumThreads(const int NumThreads)
{
    mNumThreads = NumThreads;
}

void Kernel::ImportApplication(KratosApplication::Pointer pNewApplication) {
    if (IsImported(pNewApplication->Name()))
        KRATOS_ERROR << "importing more than once the application : "
                     << pNewApplication->Name() << std::endl;

    pNewApplication->Register();
    Kernel::GetApplicationsList().insert(pNewApplication->Name());
}

std::string Kernel::Info() const { return "kernel"; }

void Kernel::PrintInfo(std::ostream& rOStream) const { rOStream << "kernel"; }

/// Print object's data.
void Kernel::PrintData(std::ostream& rOStream) const {
    rOStream << "Variables:" << std::endl;
    KratosComponents<VariableData>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Geometries:" << std::endl;
    KratosComponents<Geometry<Node<3>>>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Modelers:" << std::endl;
    KratosComponents<Modeler>().PrintData(rOStream);

    rOStream << "Loaded applications:" << std::endl;

    auto& application_list = Kernel::GetApplicationsList();
    rOStream << "number of loaded applications = " << application_list.size()
             << std::endl;
    for (auto it = application_list.begin(); it != application_list.end(); ++it)
        rOStream << "  " << *it << std::endl;
}

// To be removed with the new entry points
std::string Kernel::BuildType() {
    return GetBuildType();
}

// To be removed with the new entry points
std::string Kernel::Version() {
    return GetVersionString();
}

void Kernel::PrintParallelismSupportInfo() const
{
    #ifdef _OPENMP
    constexpr bool openmp_support = true;
    #else
    constexpr bool openmp_support = false;
    #endif

    #ifdef KRATOS_USING_MPI
    constexpr bool mpi_support = true;
    #else
    constexpr bool mpi_support = false;
    #endif

    Logger logger("");
    logger << LoggerMessage::Severity::INFO;

    if (openmp_support) {
        if (mpi_support) {
            logger << "Compiled with OpenMP and MPI support." << std::endl;
        }
        else {
            logger << "Compiled with OpenMP support." << std::endl;
        }
    }
    else if (mpi_support) {
        logger << "Compiled with MPI support." << std::endl;
    }
    else {
        logger << "Serial compilation." << std::endl;
    }

    if (openmp_support) {
        logger << "Maximum OpenMP threads: " << OpenMPUtils::GetNumThreads() << "." << std::endl;
    }

    if (mpi_support) {
        if (mIsDistributedRun) {
            const DataCommunicator& r_world = ParallelEnvironment::GetDataCommunicator("World");
            logger << "MPI world size:         " << r_world.Size() << "." << std::endl;
        }
        else {
            logger << "Running without MPI." << std::endl;
        }
    }
}

bool Kernel::mIsDistributedRun = false;

int Kernel::mNumThreads = 1;

}
