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

// Project includes
#include "internals/co_sim_io_define.h"

namespace CoSimIO {

inline void Connect(const std::string& rConnectionName, CoSimIO::SettingsType Settings);
inline void Connect(const std::string& rConnectionName, const std::string& rSettingsFileName);

inline void Disconnect(const std::string& rConnectionName);


template<class TContainerType>
inline void ImportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    int& rSize,
    TContainerType& pData);

template<class TContainerType>
inline void ExportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const int Size,
    TContainerType& pData);


template<class TDoubleContainerType, class TIntContainerType>
inline void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    int& rNumberOfNodes,
    int& rNumberOfElements,
    TDoubleContainerType& rNodalCoordinates,
    TIntContainerType& rElementConnectivities,
    TIntContainerType& rElementTypes);

template<class TDoubleContainerType, class TIntContainerType>
inline void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const int NumberOfNodes,
    const int NumberOfElements,
    TDoubleContainerType& rNodalCoordinates,
    TIntContainerType& rElementConnectivities,
    TIntContainerType& rElementTypes);


inline void IsConverged(const std::string& rConnectionName, int& rConvergenceSignal);

inline void SendControlSignal(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const CoSimIO::ControlSignal Signal);

inline void Run(const std::string& rConnectionName);

template<typename TFunctionType>
inline void Register(
    const std::string& rConnectionName,
    const std::string& rFunctionName,
    TFunctionType rFunction);

} // namespace CoSimIO

#include "internals/impl/co_sim_io_impl.h"

#endif /* KRATOS_CO_SIM_IO_H_INCLUDED */
