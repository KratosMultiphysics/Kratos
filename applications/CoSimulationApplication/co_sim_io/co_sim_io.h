// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef KRATOS_CO_SIM_IO_H_INCLUDED
#define KRATOS_CO_SIM_IO_H_INCLUDED

/*
This file defines the IO of Kratos-CoSimulation for the exchange of data with external solvers.
Only the delaration of the interface functions is defined in this file,
the corresponding definitions can be found in "impl/co_sim_io_impl.h"

By default the communication is done through files.
Support for additional means of communication can be enabled by uncommenting the following macros.
Note that this introduces dependencies such as e.g. boost (header-only version) or MPI
*/

// macros are protected to avoid redefinition in case they were defined externally (e.g. in CMake)
#ifndef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
// #define KRATOS_CO_SIM_IO_ENABLE_SOCKETS // uncomment for Sockets support
#endif
#ifndef KRATOS_CO_SIM_IO_ENABLE_MPI
// #define KRATOS_CO_SIM_IO_ENABLE_MPI // uncomment for MPI support
#endif

// System includes
#include <string>

// Project includes
#include "impl/co_sim_io_define.h"

namespace CoSimIO {

inline void Connect(const std::string& rConnectionName, CoSimIO::SettingsType Settings);
inline void Connect(const std::string& rConnectionName, const std::string& rSettingsFileName);

inline void Disconnect(const std::string& rConnectionName);


template<class TContainerType>
inline void ImportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    TContainerType& rData);

template<class TContainerType>
inline void ExportData(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    TContainerType& rData);


template<class TDoubleContainerType, class TIntContainerType>
inline void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    TDoubleContainerType& rNodalCoordinates,
    TIntContainerType& rElementConnectivities,
    TIntContainerType& rElementTypes);

template<class TDoubleContainerType, class TIntContainerType>
inline void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    TDoubleContainerType& rNodalCoordinates,
    TIntContainerType& rElementConnectivities,
    TIntContainerType& rElementTypes);


inline int IsConverged(const std::string& rConnectionName);

inline void Run(const std::string& rConnectionName);

template<typename TFunctionType>
inline void Register(
    const std::string& rConnectionName,
    const std::string& rFunctionName,
    TFunctionType rFunction);

} // namespace CoSimIO

#include "impl/co_sim_io_impl.h"

#endif // KRATOS_CO_SIM_IO_H_INCLUDED
