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

// External includes

// Project includes
#include "includes/kernel.h"
#include "includes/kratos_version.h"
#include "includes/data_communicator.h"
#include "includes/parallel_environment.h"
#include "input_output/logger.h"
#include "utilities/openmp_utils.h"

namespace Kratos {
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
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);

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

}
