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
#include "includes/registry.h"
#include "input_output/logger.h"
#include "utilities/parallel_utilities.h"

namespace Kratos {

Kernel::Kernel() 
    : mpKratosCoreApplication(Kratos::make_shared<KratosApplication>(std::string("KratosMultiphysics"))) 
{
    Initialize();
}

Kernel::Kernel(bool IsDistributedRun) 
        : mpKratosCoreApplication(Kratos::make_shared<KratosApplication>(std::string("KratosMultiphysics")))
{
    // Distributed run definition
    mIsDistributedRun = IsDistributedRun;

    // Initialize kernel
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
}

void Kernel::Initialize() {
    // Print kernel info
    this->PrintInfo();

    // Import the Kratos core application (if not already imported)
    if (!IsImported("KratosMultiphysics")) {
        // Boost is always available
        if (!Registry::HasItem("libraries.boost")) {
            Registry::AddItem<std::string>("libraries.boost");

            // When using the nonfree version of TRIANGLE, add it to the list of libraries
        #if USE_TRIANGLE_NONFREE_TPL
            Registry::AddItem<std::string>("libraries.triangle");
        #else
            // Open-source version alternative to TRIANGLE
            Registry::AddItem<std::string>("libraries.delaunator-cpp");
        #endif

            // When using the nonfree version of TETGEN, add it to the list of libraries
        #if USE_TETGEN_NONFREE_TPL
            Registry::AddItem<std::string>("libraries.tetgen");
        #endif

            // Add the libraries that are always available
            Registry::AddItem<std::string>("libraries.amgcl");
            Registry::AddItem<std::string>("libraries.benchmark");
            Registry::AddItem<std::string>("libraries.clipper");
            Registry::AddItem<std::string>("libraries.concurrentqueue");
            Registry::AddItem<std::string>("libraries.ghc");
            Registry::AddItem<std::string>("libraries.gidpost");
            Registry::AddItem<std::string>("libraries.intrusive_ptr");
            Registry::AddItem<std::string>("libraries.json");
            Registry::AddItem<std::string>("libraries.pybind11");
            Registry::AddItem<std::string>("libraries.span");
            Registry::AddItem<std::string>("libraries.tinyexpr");
            Registry::AddItem<std::string>("libraries.vexcl");
            Registry::AddItem<std::string>("libraries.zlib");
        }

        // Import core application
        this->ImportApplication(mpKratosCoreApplication);
    }
}

std::unordered_set<std::string>& Kernel::GetApplicationsList() {
    static std::unordered_set<std::string> application_list;
    return application_list;
}

std::unordered_set<std::string> Kernel::GetLibraryList() {
    std::unordered_set<std::string> library_list;

    const auto& r_item = Registry::GetItem("libraries");
    for (auto it_item = r_item.cbegin(); it_item != r_item.cend(); ++it_item) {
        library_list.insert((it_item->second)->Name());
    }

    return library_list;
}

bool Kernel::IsImported(const std::string& rApplicationName) const {
    return GetApplicationsList().find(rApplicationName) != GetApplicationsList().end();
}

bool Kernel::IsLibraryAvailable(const std::string& rLibraryName) const {
    return GetLibraryList().find(rLibraryName) != GetLibraryList().end();
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

bool Kernel::mIsDistributedRun = false;
std::string Kernel::mPyVersion = std::string("Undefined");

}
