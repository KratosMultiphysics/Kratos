//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_H_INCLUDED
#define CO_SIM_IO_H_INCLUDED

/*
This file defines the IO of Kratos-CoSimulation for the exchange of data with external solvers.
Only the delaration of the interface functions is defined in this file,
the corresponding definitions can be found in "impl/co_sim_io_impl.h"

By default the communication is done through files.
Support for additional means of communication can be enabled by uncommenting the following macros.
Note that this introduces dependencies such as e.g. boost (header-only version) or MPI
*/

// macros are protected to avoid redefinition in case they were defined externally (e.g. in CMake)
#ifndef CO_SIM_IO_USING_SOCKETS
// #define CO_SIM_IO_USING_SOCKETS // uncomment for Sockets support
#endif
#ifndef CO_SIM_IO_USING_MPI
// #define CO_SIM_IO_USING_MPI // uncomment for MPI support
#endif

// System includes
#include <string>

// Project includes
#include "impl/define.hpp"
#include "impl/info.hpp"

namespace CoSimIO {

inline Info Hello();


inline Info Connect(
    const Info& I_Settings);

inline Info Disconnect(
    const Info& I_Info);


template<class TContainerType>
inline Info ImportData(
    const Info& I_Info,
    TContainerType& rData);

template<class TContainerType>
inline Info ExportData(
    const Info& I_Info,
    const TContainerType& rData);


template<class TDoubleContainerType, class TIntContainerType>
inline Info ImportMesh(
    const Info& I_Info,
    TDoubleContainerType& rNodalCoordinates,
    TIntContainerType& rElementConnectivities,
    TIntContainerType& rElementTypes);

template<class TDoubleContainerType, class TIntContainerType>
inline Info ExportMesh(
    const Info& I_Info,
    const TDoubleContainerType& rNodalCoordinates,
    const TIntContainerType& rElementConnectivities,
    const TIntContainerType& rElementTypes);


inline Info ImportInfo(
    const Info& I_Info);

inline Info ExportInfo(
    const Info& I_Info);


inline Info IsConverged(const Info& I_Info);

inline Info Run(const Info& I_Info);

template<typename TFunctionType>
inline Info Register(
    const Info& I_Info,
    TFunctionType rFunction);

} // namespace CoSimIO

#include "impl/co_sim_io_impl.hpp"

#endif // CO_SIM_IO_H_INCLUDED
