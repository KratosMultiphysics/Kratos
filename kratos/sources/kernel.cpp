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
#include "includes/system_information.h"
#include "input_output/logger.h"
#include "utilities/parallel_utilities.h"
#include "utilities/color_utilities.h"

namespace Kratos {

Kernel::Kernel() : mpKratosCoreApplication(Kratos::make_shared<KratosApplication>(
                std::string("KratosMultiphysics"))) {
    Initialize();
}

Kernel::Kernel(bool IsDistributedRun) : mpKratosCoreApplication(Kratos::make_shared<KratosApplication>(
                std::string("KratosMultiphysics"))) {
    mIsDistributedRun = IsDistributedRun;
    Initialize();
}

Kernel::~Kernel() {
    GetApplicationsList().clear();
}

void Kernel::PrintInfo() {
    KRATOS_INFO("") << " |  /           |                  \n"
                    << " ' /   __| _` | __|  _ \\   __|    \n"
                    << " . \\  |   (   | |   (   |\\__ \\  \n"
                    << "_|\\_\\_|  \\__,_|\\__|\\___/ ____/\n"
                    << "           Multi-Physics " << Kernel::Version() << "\n"
                    << "           Compiled for "  << Kernel::OSName()  << " and " << Kernel::PythonVersion() << " with " << Kernel::Compiler() << std::endl;

    PrintParallelismSupportInfo();

    PrintSystemInfo();
}

void Kernel::Initialize() {
    this->PrintInfo();

    if (!IsImported("KratosMultiphysics")) {
        this->ImportApplication(mpKratosCoreApplication);
    }
}

std::unordered_set<std::string> &Kernel::GetApplicationsList() {
  static std::unordered_set<std::string> application_list;
  return application_list;
}

bool Kernel::IsImported(const std::string& ApplicationName) const {
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
    KratosComponents<Geometry<Node>>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Modelers:" << std::endl;
    KratosComponents<Modeler>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Loaded applications:" << std::endl;
    auto& application_list = Kernel::GetApplicationsList();
    rOStream << "    Number of loaded applications = " << application_list.size() << std::endl;
    for (auto it = application_list.begin(); it != application_list.end(); ++it)
        rOStream << "    " << *it << std::endl;
}

// To be removed with the new entry points
std::string Kernel::BuildType() {
    return GetBuildType();
}

// To be removed with the new entry points
std::string Kernel::Version() {
    return GetVersionString();
}

std::string Kernel::OSName() {
    return GetOSName();
}

std::string Kernel::PythonVersion() {
    return mPyVersion;
}

std::string Kernel::Compiler() {
    return GetCompiler();
}

void Kernel::SetPythonVersion(std::string pyVersion) {
    mPyVersion = pyVersion;
}

void Kernel::PrintParallelismSupportInfo() const
{
    #ifdef KRATOS_SMP_NONE
    constexpr bool threading_support = false;
    #else
    constexpr bool threading_support = true;
    #endif

    #ifdef KRATOS_USING_MPI
    constexpr bool mpi_support = true;
    #else
    constexpr bool mpi_support = false;
    #endif

    Logger logger("");
    logger << LoggerMessage::Severity::INFO;

    if (threading_support) {
        if (mpi_support) {
            logger << "Compiled with threading and MPI support." << std::endl;
        }
        else {
            logger << "Compiled with threading support." << std::endl;
        }
    }
    else if (mpi_support) {
        logger << "Compiled with MPI support." << std::endl;
    }
    else {
        logger << "Serial compilation." << std::endl;
    }

    if (threading_support) {
        logger << "Maximum number of threads: " << ParallelUtilities::GetNumThreads() << "." << std::endl;
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

void Kernel::PrintSystemInfo() const
{
    // TODO: Implement these methods
    // const std::string hyper_threading = SystemInformation::CpuHyperThreading() ? "Enabled" : "Disabled";
    // KRATOS_INFO("") << BOLDFONT(FWHT("System information:"))                                                                                                         << RST << "\n"
    //                 << ITAFONT(FWHT("\tOS version:\t\t"))        << KCYN << SystemInformation::OSVersion()                                                           << RST << "\n"
    //                 << ITAFONT(FWHT("\tCPU architecture:\t"))    << KGRN << SystemInformation::CpuArchitecture()                                                     << RST << "\n"
    //                 << ITAFONT(FWHT("\tCPU logical cores:\t"))   << KGRN << std::to_string(SystemInformation::CpuLogicalCores())                                     << RST << "\n"
    //                 << ITAFONT(FWHT("\tCPU physical cores:\t"))  << KGRN << std::to_string(SystemInformation::CpuPhysicalCores())                                    << RST << "\n"
    //                 << ITAFONT(FWHT("\tCPU clock speed:\t"))     << KGRN << SystemInformation::GenerateClockSpeed(SystemInformation::CpuClockSpeed())                << RST << "\n"
    //                 << ITAFONT(FWHT("\tCPU Hyper-Threading:\t")) << KGRN << hyper_threading                                                                          << RST << "\n"
    //                 << ITAFONT(FWHT("\tRAM total:\t\t"))         << KYEL << SystemInformation::GenerateDataSize(SystemInformation::RamTotal())                       << RST << "\n"
    //                 << ITAFONT(FWHT("\tRAM free:\t\t"))          << KYEL << SystemInformation::GenerateDataSize(SystemInformation::RamFree())                        << RST << std::endl;
}

bool Kernel::mIsDistributedRun = false;
std::string Kernel::mPyVersion = std::string("Undefined");

}
