//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "shock_detection_process.h"

namespace Kratos
{

/* Public functions *******************************************************/

void ShockDetectionProcess::ExecuteInitialize()
{
    // Calculate the nodal area
    CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>(
        mrModelPart,
        mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE)).Execute();

    // Calculate the nodal neighbours for the edge-based shock capturing
    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
    FindGlobalNodalNeighboursProcess(r_data_communicator, mrModelPart).Execute();
}

void ShockDetectionProcess::EdgeBasedShockDetection(
    const Variable<double>& rShockVariable,
    const Variable<array_1d<double, 3>>& rShockGradientVariable)
{
    // Specialize the edge based shock detection
    EdgeBasedShockDetectionSpecialization<>(rShockVariable, rShockGradientVariable);
}

void ShockDetectionProcess::EdgeBasedShockDetection(
    const VariableComponentType& rShockVariable,
    const Variable<array_1d<double, 3>>& rShockGradientVariable)
{
    // Specialize the edge based shock detection
    EdgeBasedShockDetectionSpecialization<>(rShockVariable, rShockGradientVariable);
}

std::string ShockDetectionProcess::Info() const
{
    std::stringstream buffer;
    buffer << "ShockDetectionProcess";
    return buffer.str();
}

void ShockDetectionProcess::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "ShockDetectionProcess";
}

void ShockDetectionProcess::PrintData(std::ostream &rOStream) const
{
    mrModelPart.PrintData(rOStream);
    rOStream << "Update nodal area at each step: " << mUpdateNodalAreaAtEachStep << std::endl;
    rOStream << "Update nodal neighbours at each step: " << mUpdateNodalNeighboursAtEachStep << std::endl;
    rOStream << "Shock sensor variable: " << mrShockSensorVariable.Name() << std::endl;
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/

};  // namespace Kratos.
