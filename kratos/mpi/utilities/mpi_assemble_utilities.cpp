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
    const TMap<int>& rNodalValuesMap) const
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
    const TMap<double>& rNodalValuesMap) const
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
    const TMap<array_1d<double, 3>>& rNodalValuesMap) const
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
    const TMap<int>& rNodalValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateNonHistoricalValue<int, NodeType>);
    r_communicator.SynchronizeNonHistoricalVariable(rVariable);
}

void MPIAssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rNodalValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateNonHistoricalValue<double, NodeType>);
    r_communicator.SynchronizeNonHistoricalVariable(rVariable);
}

void MPIAssembleUtilities::AssembleNonHistoricalNodalDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rNodalValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Nodes();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rNodalValuesMap,
        BaseType::UpdateNonHistoricalValue<array_1d<double, 3>, NodeType>);
    r_communicator.SynchronizeNonHistoricalVariable(rVariable);
}

void MPIAssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const TMap<int>& rElementValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Elements();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rElementValuesMap,
        BaseType::UpdateNonHistoricalValue<int, ElementType>);
}

void MPIAssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const TMap<double>& rElementValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Elements();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rElementValuesMap,
        BaseType::UpdateNonHistoricalValue<double, ElementType>);
}

void MPIAssembleUtilities::AssembleElementDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rElementValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Elements();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rElementValuesMap,
        BaseType::UpdateNonHistoricalValue<array_1d<double, 3>, ElementType>);
}

void MPIAssembleUtilities::AssembleConditionDataWithValuesMap(ModelPart& rModelPart,
                                                           const Variable<int>& rVariable,
                                                           const TMap<int>& rConditionValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Conditions();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rConditionValuesMap,
        BaseType::UpdateNonHistoricalValue<int, ConditionType>);
}

void MPIAssembleUtilities::AssembleConditionDataWithValuesMap(ModelPart& rModelPart,
                                                           const Variable<double>& rVariable,
                                                           const TMap<double>& rConditionValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Conditions();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rConditionValuesMap,
        BaseType::UpdateNonHistoricalValue<double, ConditionType>);
}

void MPIAssembleUtilities::AssembleConditionDataWithValuesMap(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const TMap<array_1d<double, 3>>& rConditionValuesMap) const
{
    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_container = r_communicator.LocalMesh().Conditions();
    MPIAssembleUtilities::MPIAssembleDataWithEntityValuesMap(
        r_container, r_communicator.GetDataCommunicator(), rVariable, rConditionValuesMap,
        BaseType::UpdateNonHistoricalValue<array_1d<double, 3>, ConditionType>);
}

} // namespace Kratos