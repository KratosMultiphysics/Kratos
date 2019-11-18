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

namespace CoSimIO {

static void Connect(const char* pConnectionName, const char* pSettingsFileName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF(HasIO(pConnectionName)) << "A connection for " << pConnectionName << " already exists!" << std::endl;

    s_co_sim_ios[std::string(pConnectionName)] = std::unique_ptr<CoSimIOImpl>(new CoSimIOImpl(pConnectionName, pSettingsFileName)); // make_unique is C++14
    GetIO(pConnectionName).Connect();
}

static void Disconnect(const char* pConnectionName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(pConnectionName)) << "Trying to disconnect connection " << pConnectionName << " which does not exist!" << std::endl;

    GetIO(pConnectionName).Disconnect();
    s_co_sim_ios.erase(std::string(pConnectionName));
}

static bool IsConverged(const char* pConnectionName)
{
    return Internals::GetIO(pConnectionName).IsConverged();
}

static void ImportData(
    const char* pConnectionName,
    const char* pIdentifier,
    int* pSize,
    double** ppData)
{
    Internals::GetIO(pConnectionName).ImportData(pIdentifier, pSize, ppData);
}

static void ExportData(
    const char* pConnectionName,
    const char* pIdentifier,
    const int Size,
    const double* pData)
{
    Internals::GetIO(pConnectionName).ExportData(pIdentifier, Size, pData);
}

static void ImportMesh(
    const char* pConnectionName,
    const char* pIdentifier,
    int* pNumberOfNodes,
    int* pNumberOfElements,
    double** ppNodalCoordinates,
    int** ppElementConnectivities,
    int** ppElementTypes)
{
    Internals::GetIO(pConnectionName).ImportMesh(pIdentifier, pNumberOfNodes, pNumberOfElements, ppNodalCoordinates, ppElementConnectivities, ppElementTypes);
}

static void ExportMesh(
    const char* pConnectionName,
    const char* pIdentifier,
    const int NumberOfNodes,
    const int NumberOfElements,
    const double* pNodalCoordinates,
    const int* pElementConnectivities,
    const int* pElementTypes)
{
    Internals::GetIO(pConnectionName).ExportMesh(pIdentifier, NumberOfNodes, NumberOfElements, pNodalCoordinates, pElementConnectivities, pElementTypes);
}

static void ImportGeometry(const char* pConnectionName)
{
    Internals::GetIO(pConnectionName).ImportGeometry();
}

static void ExportGeometry(const char* pConnectionName)
{
    Internals::GetIO(pConnectionName).ExportGeometry();
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

} // namespace CoSimIO


#endif /* KRATOS_CO_SIM_IO_H_INCLUDED */
