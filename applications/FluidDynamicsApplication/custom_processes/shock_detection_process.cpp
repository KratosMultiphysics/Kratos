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
    // Initialize the non-historical database variables
    // Note that this very first initialization hast to be done out of a parallel region
    for (auto& r_node : mrModelPart.Nodes()) {
        r_node.SetValue(NEIGHBOUR_NODES, NEIGHBOUR_NODES.Zero());
        r_node.SetValue((*mpShockSensorVariable), (*mpShockSensorVariable).Zero());
    }
}

void ShockDetectionProcess::ExecuteInitializeSolutionStep()
{
    // Calculate the nodal area
    if (!mNodalAreaAlreadyComputed || mUpdateNodalAreaAtEachStep) {
        CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>(
            mrModelPart,
            mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE)).Execute();
        mNodalAreaAlreadyComputed = true;
    }

    // Calculate the nodal neighbours for the edge-based shock capturing
    if (!mNodalNeighboursAlreadyComputed || mUpdateNodalNeighboursAtEachStep) {
        const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
        FindGlobalNodalNeighboursProcess(r_data_communicator, mrModelPart).Execute();
        mNodalNeighboursAlreadyComputed = true;
    }

    // Perform the edge based shock detection
    if (mShockVariableIsDouble) {
        EdgeBasedShockDetectionSpecialization<>(*mpShockDoubleVariable, *mpShockGradientVariable);
    } else {
        EdgeBasedShockDetectionSpecialization<>(*mpShockComponentVariable, *mpShockGradientVariable);
    }
}

void ShockDetectionProcess::Execute()
{
    // Initialize the shock sensor, nodal mass and neighbours
    ExecuteInitialize();

    // Perform the edge based shock detection
    ExecuteInitializeSolutionStep();
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
    const auto& r_shock_sensor_variable = KratosComponents<Variable<double>>::Get(mpShockSensorVariable->Name());
    rOStream << "Shock sensor variable: " << r_shock_sensor_variable.Name() << std::endl;
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/

};  // namespace Kratos.
