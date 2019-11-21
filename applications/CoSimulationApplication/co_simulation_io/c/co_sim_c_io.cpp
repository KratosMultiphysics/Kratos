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

// Project includes
extern "C" {
#include "co_sim_c_io.h"
}
#include "../co_sim_io.h"


void CoSimIO_Connect(const char* pConnectionName, const char* pSettingsFileName)
{
    CoSimIO::Connect(pConnectionName, pSettingsFileName);
}

void CoSimIO_Disconnect(const char* pConnectionName)
{
    CoSimIO::Disconnect(pConnectionName);
}

void CoSimIO_ImportData(
    const char* pConnectionName,
    const char* pIdentifier,
    int* pSize,
    double** ppData)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerRawMemory<double>(ppData, *pSize)); // make_unique is c++14
    CoSimIO::ImportData(pConnectionName, pIdentifier, p_container.get());
    // TODO check return of Size!
}

void CoSimIO_ExportData(
    const char* pConnectionName,
    const char* pIdentifier,
    int Size,
    double* pData)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerRawMemory<double>(&pData, Size)); // make_unique is c++14
    CoSimIO::ExportData(pConnectionName, pIdentifier, p_container.get());
}

void ImportMesh(
    const char* pConnectionName,
    const char* pIdentifier,
    int* pNumberOfNodes,
    int* pNumberOfElements,
    double** ppNodalCoordinates,
    int** ppElementConnectivities,
    int** ppElementTypes)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container_coords(new DataContainerRawMemory<double>(ppNodalCoordinates, *pNumberOfNodes)); // make_unique is c++14
    std::unique_ptr<DataContainer<int>> p_container_conn(new DataContainerRawMemory<int>(ppElementConnectivities, *pNumberOfElements)); // make_unique is c++14
    std::unique_ptr<DataContainer<int>> p_container_types(new DataContainerRawMemory<int>(ppElementTypes, *pNumberOfElements)); // make_unique is c++14
    CoSimIO::ImportMesh(pConnectionName, pIdentifier, p_container_coords.get(), p_container_conn.get(), p_container_types.get());
    // TODO check return of Sizes!
}

void CoSimIO_ExportMesh(
    const char* pConnectionName,
    const char* pIdentifier,
    int NumberOfNodes,
    int NumberOfElements,
    double* pNodalCoordinates,
    int* pElementConnectivities,
    int* pElementTypes)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container_coords(new DataContainerRawMemory<double>(&pNodalCoordinates, NumberOfNodes)); // make_unique is c++14
    std::unique_ptr<DataContainer<int>> p_container_conn(new DataContainerRawMemory<int>(&pElementConnectivities, NumberOfElements)); // make_unique is c++14
    std::unique_ptr<DataContainer<int>> p_container_types(new DataContainerRawMemory<int>(&pElementTypes, NumberOfElements)); // make_unique is c++14
    CoSimIO::ExportMesh(pConnectionName, pIdentifier, p_container_coords.get(), p_container_conn.get(), p_container_types.get());
}
