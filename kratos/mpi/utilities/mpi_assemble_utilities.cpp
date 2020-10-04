//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <unordered_map>

// External includes

// Project includes

// Include base h
#include "mpi/utilities/mpi_assemble_utilities.h"

namespace Kratos
{
void MPIAssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<NodeType, int>& rNodalValuesMap) const
{
    BaseType::CheckHistoricalVariable(rModelPart, rVariable);

    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateHistoricalNodalValue<int>);
    r_communicator.SynchronizeVariable(rVariable);
}

void MPIAssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<NodeType, double>& rNodalValuesMap) const
{
    BaseType::CheckHistoricalVariable(rModelPart, rVariable);

    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateHistoricalNodalValue<double>);
    r_communicator.SynchronizeVariable(rVariable);
}

void MPIAssembleUtilities::AssembleCurrentDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const
{
    BaseType::CheckHistoricalVariable(rModelPart, rVariable);

    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateHistoricalNodalValue<array_1d<double, 3>>);
    r_communicator.SynchronizeVariable(rVariable);
}

void MPIAssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<NodeType, int>& rNodalValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateNonHistoricalValue<NodeType, int>);
    r_communicator.SynchronizeNonHistoricalVariable(rVariable);
}

void MPIAssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<NodeType, double>& rNodalValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateNonHistoricalValue<NodeType, double>);
    r_communicator.SynchronizeNonHistoricalVariable(rVariable);
}

void MPIAssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateNonHistoricalValue<NodeType, array_1d<double, 3>>);
    r_communicator.SynchronizeNonHistoricalVariable(rVariable);
}

void MPIAssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TGPMap<ElementType, int>& rElementValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Elements();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rElementValuesMap,
        BaseType::UpdateNonHistoricalValue<ElementType, int>);
}

void MPIAssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TGPMap<ElementType, double>& rElementValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Elements();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rElementValuesMap,
        BaseType::UpdateNonHistoricalValue<ElementType, double>);
}

void MPIAssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<ElementType, array_1d<double, 3>>& rElementValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Elements();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rElementValuesMap,
        BaseType::UpdateNonHistoricalValue<ElementType, array_1d<double, 3>>);
}

void MPIAssembleUtilities::AssembleConditionDataWithValuesMap(ModelPart& rModelPart,
                                                           const Variable<int>& rVariable,
                                                           const TGPMap<ConditionType, int>& rConditionValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Conditions();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rConditionValuesMap,
        BaseType::UpdateNonHistoricalValue<ConditionType, int>);
}

void MPIAssembleUtilities::AssembleConditionDataWithValuesMap(ModelPart& rModelPart,
                                                           const Variable<double>& rVariable,
                                                           const TGPMap<ConditionType, double>& rConditionValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Conditions();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rConditionValuesMap,
        BaseType::UpdateNonHistoricalValue<ConditionType, double>);
}

void MPIAssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TGPMap<ConditionType, array_1d<double, 3>>& rConditionValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Conditions();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rConditionValuesMap,
        BaseType::UpdateNonHistoricalValue<ConditionType, array_1d<double, 3>>);
}

} // namespace Kratos