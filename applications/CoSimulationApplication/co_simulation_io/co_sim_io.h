// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_IO_H_INCLUDED
#define KRATOS_CO_SIM_IO_H_INCLUDED

/*
This file defines the IO of Kratos-CoSimulation for the exchange of data
with external solvers.
By default the communication is done through files,
support for sockets and MPI can optionally be enabled
*/

// #define KRATOS_CO_SIM_IO_ENABLE_SOCKETS // uncomment for Sockets support
// #define KRATOS_CO_SIM_IO_ENABLE_MPI // uncomment for MPI support

// System includes
#include <string>
#include <memory>

// Project includes
#include "impl/co_sim_io_impl.h"

namespace CoSim {

using Internals::CoSimIO; // TODO remove this

static void Connect(const char* pName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF(HasIO(pName)) << "A CoSimIO for " << pName << " already exists!" << std::endl;

    s_co_sim_ios[std::string(pName)] = std::unique_ptr<CoSimIO>(new CoSimIO("rName", "rSettings")); // make_unique is C++14
    GetIO(pName).Connect();
}

static void Disconnect(const char* pName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(pName)) << "Trying to disconnect CoSimIO " << pName << " which does not exist!" << std::endl;

    GetIO(pName).Disconnect();
    s_co_sim_ios.erase(std::string(pName));
}

static bool IsConverged(const char* pName)
{
    return Internals::GetIO(pName).IsConverged();
}

static void ImportData(
    const char* pName,
    const char* pIdentifier,
    int* pSize,
    double** ppData)
{
    Internals::GetIO(pName).ImportData(pIdentifier, pSize, ppData);
}

static void ExportData(
    const char* pName,
    const char* pIdentifier,
    const int Size,
    const double* pData)
{
    Internals::GetIO(pName).ExportData(pIdentifier, Size, pData);
}

static void ImportMesh(
    const char* pName,
    const char* pIdentifier,
    int* pNumberOfNodes,
    int* pNumberOfElements,
    double** ppNodalCoordinates,
    int** ppElementConnectivities,
    int** ppElementTypes)
{
    Internals::GetIO(pName).ImportMesh(pIdentifier, pNumberOfNodes, pNumberOfElements, ppNodalCoordinates, ppElementConnectivities, ppElementTypes);
}

static void ExportMesh(
    const char* pName,
    const char* pIdentifier,
    const int NumberOfNodes,
    const int NumberOfElements,
    const double* pNodalCoordinates,
    const int* pElementConnectivities,
    const int* pElementTypes)
{
    Internals::GetIO(pName).ExportMesh(pIdentifier, NumberOfNodes, NumberOfElements, pNodalCoordinates, pElementConnectivities, pElementTypes);
}

static void ImportGeometry(const char* pName)
{
    Internals::GetIO(pName).ImportGeometry();
}

static void ExportGeometry(const char* pName)
{
    Internals::GetIO(pName).ExportGeometry();
}

namespace {
    // defining aliases to the functions to avoid the "unused-functions" compiler-warning if not all functions are used
    const auto _alias_Connect = Connect;
    const auto _alias_Disconnect = Disconnect;
    const auto _alias_IsConverged = IsConverged;
    const auto _alias_ImportData = ImportData;
    const auto _alias_ExportData = ExportData;
    const auto _alias_ImportMesh = ImportMesh;
    const auto _alias_ExportMesh = ExportMesh;
    const auto _alias_ImportGeometry = ImportGeometry;
    const auto _alias_ExportGeometry = ExportGeometry;
}

} // namespace CoSim


#endif /* KRATOS_CO_SIM_IO_H_INCLUDED */
