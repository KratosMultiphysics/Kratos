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
#include "co_sim_io_impl_2.h"

namespace CoSimIO {

namespace Internals {
// TODO make sure this is unique even across compilation units (test somehow)
static std::unordered_map<std::string, std::unique_ptr<CoSimIOImpl>> s_co_sim_ios;

static bool HasIO(const std::string& rConnectionName)
{
    return s_co_sim_ios.find(rConnectionName) != s_co_sim_ios.end();
}

static CoSimIOImpl& GetIO(const std::string& rConnectionName)
{
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(rConnectionName)) << "Trying to use connection " << rConnectionName << " which does not exist!" << std::endl;
    return *s_co_sim_ios.at(rConnectionName);
}

} // namespace Internals

inline void Connect(const std::string& rConnectionName, const std::string& pSettingsFileName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF(HasIO(rConnectionName)) << "A connection for " << rConnectionName << " already exists!" << std::endl;

    s_co_sim_ios[std::string(rConnectionName)] = std::unique_ptr<CoSimIOImpl>(new CoSimIOImpl(rConnectionName, pSettingsFileName));
    GetIO(rConnectionName).Connect();
}

inline void Disconnect(const std::string& rConnectionName)
{
    using namespace Internals;
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(rConnectionName)) << "Trying to disconnect connection " << rConnectionName << " which does not exist!" << std::endl;

    GetIO(rConnectionName).Disconnect();
    s_co_sim_ios.erase(std::string(rConnectionName));
}

// Version for C++, there this input is a std::vector, which we have to wrap before passing it on
template<>
inline void ImportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    std::vector<double>& pData)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerStdVector<double>(pData));
    GetIO(rConnectionName).ImportData(rIdentifier, *p_container);
}

// Version for C and fortran, there we already get a container
template<>
inline void ImportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    CoSimIO::Internals::DataContainer<double>& pData)
{
    Internals::GetIO(rConnectionName).ImportData(rIdentifier, pData);
}

// Version for C++, there this input is a std::vector, which we have to wrap before passing it on
template<>
inline void ExportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    std::vector<double>& pData)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container(new DataContainerStdVector<double>(pData));
    GetIO(rConnectionName).ExportData(rIdentifier, *p_container);
}

// Version for C and fortran, there we already get a container
template<>
inline void ExportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    CoSimIO::Internals::DataContainer<double>& pData)
{
    Internals::GetIO(rConnectionName).ExportData(rIdentifier, pData);
}

template<>
inline void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    std::vector<double>& pNodalCoordinates,
    std::vector<int>& pElementConnectivities,
    std::vector<int>& pElementTypes)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container_coords(new DataContainerStdVector<double>(pNodalCoordinates));
    std::unique_ptr<DataContainer<int>> p_container_conn(new DataContainerStdVector<int>(pElementConnectivities));
    std::unique_ptr<DataContainer<int>> p_container_types(new DataContainerStdVector<int>(pElementTypes));
    Internals::GetIO(rConnectionName).ImportMesh(rIdentifier, *p_container_coords, *p_container_conn, *p_container_types);
}

template<>
inline void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    CoSimIO::Internals::DataContainer<double>& pNodalCoordinates,
    CoSimIO::Internals::DataContainer<int>& pElementConnectivities,
    CoSimIO::Internals::DataContainer<int>& pElementTypes)
{
    Internals::GetIO(rConnectionName).ImportMesh(rIdentifier, pNodalCoordinates, pElementConnectivities, pElementTypes);
}

template<>
inline void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    std::vector<double>& pNodalCoordinates,
    std::vector<int>& pElementConnectivities,
    std::vector<int>& pElementTypes)
{
    using namespace CoSimIO::Internals;
    std::unique_ptr<DataContainer<double>> p_container_coords(new DataContainerStdVector<double>(pNodalCoordinates));
    std::unique_ptr<DataContainer<int>> p_container_conn(new DataContainerStdVector<int>(pElementConnectivities));
    std::unique_ptr<DataContainer<int>> p_container_types(new DataContainerStdVector<int>(pElementTypes));
    Internals::GetIO(rConnectionName).ExportMesh(rIdentifier, *p_container_coords, *p_container_conn, *p_container_types);
}

template<>
inline void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    CoSimIO::Internals::DataContainer<double>& pNodalCoordinates,
    CoSimIO::Internals::DataContainer<int>& pElementConnectivities,
    CoSimIO::Internals::DataContainer<int>& pElementTypes)
{
    Internals::GetIO(rConnectionName).ExportMesh(rIdentifier, pNodalCoordinates, pElementConnectivities, pElementTypes);
}

inline bool IsConverged(const std::string& rConnectionName)
{
    return Internals::GetIO(rConnectionName).IsConverged();
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
    Internals::GetIO(rConnectionName).Register(rFunctionName, rFunction);
}

} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_IO_IMPL_H_INCLUDED */
