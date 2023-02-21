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

#ifndef CO_SIM_IO_INCLUDED
#define CO_SIM_IO_INCLUDED

/*
This file defines the IO of Kratos-CoSimulation for the exchange of data with external solvers.
*/

// System includes

// Project includes
#include "includes/define.hpp"
#include "includes/info.hpp"
#include "includes/model_part.hpp"
#include "includes/version.hpp"

namespace CoSimIO {

Info CO_SIM_IO_API Hello();


Info CO_SIM_IO_API Connect(
    const Info& I_Settings);

Info CO_SIM_IO_API Disconnect(
    const Info& I_Info);


template<class TContainerType>
Info CO_SIM_IO_API ImportData(
    const Info& I_Info,
    TContainerType& rData);

template<class TContainerType>
Info CO_SIM_IO_API ExportData(
    const Info& I_Info,
    const TContainerType& rData);


Info CO_SIM_IO_API ImportMesh(
    const Info& I_Info,
    ModelPart& O_ModelPart);

Info CO_SIM_IO_API ExportMesh(
    const Info& I_Info,
    const ModelPart& I_ModelPart);

Info CO_SIM_IO_API ImportInfo(
    const Info& I_Info);

Info CO_SIM_IO_API ExportInfo(
    const Info& I_Info);


Info CO_SIM_IO_API Run(const Info& I_Info);

template<typename TFunctionType>
Info CO_SIM_IO_API Register(
    const Info& I_Info,
    TFunctionType rFunction);

} // namespace CoSimIO

#endif // CO_SIM_IO_INCLUDED
