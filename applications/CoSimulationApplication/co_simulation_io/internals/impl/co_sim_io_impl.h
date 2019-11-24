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

#ifndef KRATOS_CO_SIM_IO_IMPL_H_INCLUDED
#define KRATOS_CO_SIM_IO_IMPL_H_INCLUDED

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
#include "../co_sim_connection.h"

namespace CoSimIO {

namespace Internals {
// TODO make sure this is unique even across compilation units (test somehow)
static std::unordered_map<std::string, std::unique_ptr<CoSimConnection>> s_co_sim_ios;

static bool HasIO(const std::string& rConnectionName)
{
    return s_co_sim_ios.find(rConnectionName) != s_co_sim_ios.end();
}

static CoSimConnection& GetIO(const std::string& rConnectionName)
{
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(rConnectionName)) << "Trying to use connection \"" << rConnectionName << "\" which does not exist!" << std::endl;
    return *s_co_sim_ios.at(rConnectionName);
}

} // namespace Internals

inline void Connect(const std::string& rConnectionName, const std::string& pSettingsFileName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF(HasIO(rConnectionName)) << "A connection for \"" << rConnectionName << "\" already exists!" << std::endl;

    s_co_sim_ios[std::string(rConnectionName)] = std::unique_ptr<CoSimConnection>(new CoSimConnection(rConnectionName, pSettingsFileName));
    GetIO(rConnectionName).Connect();
}

inline void Disconnect(const std::string& rConnectionName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(rConnectionName)) << "Trying to disconnect connection \"" << rConnectionName << "\" which does not exist!" << std::endl;

    GetIO(rConnectionName).Disconnect();
    s_co_sim_ios.erase(std::string(rConnectionName));
}

// Version for C++, there this input is a std::vector, which we have to wrap before passing it on
template<>
inline void ImportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    int& rSize,
    std::vector<double>& pData)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerStdVector<double>(pData));
    GetIO(rConnectionName).ImportData(rIdentifier, rSize, *p_container);
}

// Version for C and fortran, there we already get a container
template<>
inline void ImportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    int& rSize,
    CoSimIO::Internals::DataContainer<double>& pData)
{
    Internals::GetIO(rConnectionName).ImportData(rIdentifier, rSize, pData);
}

// Version for C++, there this input is a std::vector, which we have to wrap before passing it on
template<>
inline void ExportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const int Size,
    std::vector<double>& pData)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerStdVector<double>(pData));
    GetIO(rConnectionName).ExportData(rIdentifier, Size, *p_container);
}

// Version for C and fortran, there we already get a container
template<>
inline void ExportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const int Size,
    CoSimIO::Internals::DataContainer<double>& pData)
{
    Internals::GetIO(rConnectionName).ExportData(rIdentifier, Size, pData);
}

template<>
inline void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    int& rNumberOfNodes,
    int& rNumberOfElements,
    std::vector<double>& pNodalCoordinates,
    std::vector<int>& pElementConnectivities,
    std::vector<int>& pElementTypes)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container_coords(new DataContainerStdVector<double>(pNodalCoordinates));
    std::unique_ptr<DataContainer<int>> p_container_conn(new DataContainerStdVector<int>(pElementConnectivities));
    std::unique_ptr<DataContainer<int>> p_container_types(new DataContainerStdVector<int>(pElementTypes));
    Internals::GetIO(rConnectionName).ImportMesh(rIdentifier, rNumberOfNodes, rNumberOfElements, *p_container_coords, *p_container_conn, *p_container_types);
}

template<>
inline void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    int& rNumberOfNodes,
    int& rNumberOfElements,
    CoSimIO::Internals::DataContainer<double>& pNodalCoordinates,
    CoSimIO::Internals::DataContainer<int>& pElementConnectivities,
    CoSimIO::Internals::DataContainer<int>& pElementTypes)
{
    Internals::GetIO(rConnectionName).ImportMesh(rIdentifier, rNumberOfNodes, rNumberOfElements, pNodalCoordinates, pElementConnectivities, pElementTypes);
}

template<>
inline void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const int NumberOfNodes,
    const int NumberOfElements,
    std::vector<double>& pNodalCoordinates,
    std::vector<int>& pElementConnectivities,
    std::vector<int>& pElementTypes)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container_coords(new DataContainerStdVector<double>(pNodalCoordinates));
    std::unique_ptr<DataContainer<int>> p_container_conn(new DataContainerStdVector<int>(pElementConnectivities));
    std::unique_ptr<DataContainer<int>> p_container_types(new DataContainerStdVector<int>(pElementTypes));
    Internals::GetIO(rConnectionName).ExportMesh(rIdentifier, NumberOfNodes, NumberOfElements, *p_container_coords, *p_container_conn, *p_container_types);
}

template<>
inline void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const int NumberOfNodes,
    const int NumberOfElements,
    CoSimIO::Internals::DataContainer<double>& pNodalCoordinates,
    CoSimIO::Internals::DataContainer<int>& pElementConnectivities,
    CoSimIO::Internals::DataContainer<int>& pElementTypes)
{
    Internals::GetIO(rConnectionName).ExportMesh(rIdentifier, NumberOfNodes, NumberOfElements, pNodalCoordinates, pElementConnectivities, pElementTypes);
}

inline void IsConverged(const std::string& rConnectionName, int* pConvergenceSignal)
{
    Internals::GetIO(rConnectionName).IsConverged(pConvergenceSignal);
}

inline void Run(const std::string& rConnectionName)
{
    Internals::GetIO(rConnectionName).Run();
}

template<typename TFunctionType>
inline void Register(
    const std::string& rConnectionName,
    const std::string& rFunctionName,
    TFunctionType rFunction)
{
    // TODO If I use structs (base-class etc, constructed right here) then I can save all of them in one map and do the corresponding actions with them, since I know thich interface-fct to call!
    // => this means that I have to specialize this fct, but the CoSimConnection class will be much cleaner
    // (maybe internally they could be registered as void* to void* but probably not (a good idea))
    Internals::GetIO(rConnectionName).Register(rFunctionName, rFunction);
}

} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_IO_IMPL_H_INCLUDED */
