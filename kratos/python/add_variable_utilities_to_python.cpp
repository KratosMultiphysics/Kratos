//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "python/add_variable_utilities_to_python.h"
#include "includes/define_python.h"
#include "processes/process.h"

// Variable utilities
#include "utilities/variable_utils.h"


namespace Kratos {
namespace Python {

void VariableUtilsUpdateCurrentPosition(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes);
}

void VariableUtilsUpdateCurrentPositionWithVariable(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes,
    const VariableUtils::ArrayVarType &rUpdateVariable
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes, rUpdateVariable);
}

void VariableUtilsUpdateCurrentPositionWithVariableAndPosition(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes,
    const VariableUtils::ArrayVarType &rUpdateVariable,
    const IndexType BufferPosition
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes, rUpdateVariable, BufferPosition);
}

template<class TVarType>
void VariableUtilsCopyModelPartNodalVar(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void VariableUtilsCopyModelPartNodalVarWithDestination(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TVarType &rDestinationVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void CopyModelPartNodalVarToNonHistoricalVar(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void CopyModelPartNodalVarToNonHistoricalVarWithDestination(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TVarType &rDestinationVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}



void AddVariableUtilitiesToPython(pybind11::module &m)  
{

    namespace py = pybind11; 

    py::class_<VariableUtils>(m, "VariableUtils")
        .def(py::init<>())
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<bool>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<double>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double,4>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double,6>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double,9>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Matrix>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<bool>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<double>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double,3>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double,4>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double,6>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double,9>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Vector>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Matrix>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<bool>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<double>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double,4>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double,6>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double,9>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Matrix>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<bool>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<double>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double,3>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double,4>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double,6>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double,9>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Vector>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Matrix>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<bool>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<double>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,4>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,6>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,9>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Vector>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Matrix>>)
        .def("SetVectorVar", &VariableUtils::SetVectorVar)
        .def("SetVectorVar", &VariableUtils::SetVectorVarForFlag)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<Variable<double>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<Variable<double>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetNonHistoricalVectorVar", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetVariable", &VariableUtils::SetVariable<int>)
        .def("SetVariable", &VariableUtils::SetVariable<bool>)
        .def("SetVariable", &VariableUtils::SetVariable<double>)
        .def("SetVariable", &VariableUtils::SetVariable<array_1d<double, 3>>)
        .def("SetVariable", &VariableUtils::SetVariable<Vector>)
        .def("SetVariable", &VariableUtils::SetVariable<Matrix>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<bool>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<double>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 3>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 4>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 6>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 9>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<Vector>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<Matrix>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<int>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<double>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<array_1d<double, 3>>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Matrix, ModelPart::ElementsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::NodesContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ConditionsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::NodesContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ConditionsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::NodesContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::ConditionsContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::ElementsContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::NodesContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::ConditionsContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::ElementsContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("SaveVectorVar", &VariableUtils::SaveVectorVar)
        .def("SaveScalarVar", &VariableUtils::SaveScalarVar)
        .def("SaveVectorNonHistoricalVar", &VariableUtils::SaveVectorNonHistoricalVar)
        .def("SaveScalarNonHistoricalVar", &VariableUtils::SaveScalarNonHistoricalVar)
        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        .def("CopyVectorVar", &VariableUtils::CopyVectorVar)
        .def("CopyComponentVar", &VariableUtils::CopyComponentVar)
        .def("CopyScalarVar", &VariableUtils::CopyScalarVar)
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<double> >)
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 4> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 6> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 9> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 3> > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 4> > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 6> > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 9> > > )
        .def("ApplyFixity", &VariableUtils::ApplyFixity<Variable<double>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<Variable<double>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalNodeVectorVariable)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumNonHistoricalNodeVectorVariable", &VariableUtils::SumNonHistoricalNodeVectorVariable)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<Variable<double>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumConditionVectorVariable", &VariableUtils::SumConditionVectorVariable)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<Variable<double>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumElementVectorVariable", &VariableUtils::SumElementVectorVariable)
        .def("AddDof", &VariableUtils::AddDof<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("CheckVariableKeys", &VariableUtils::CheckVariableKeys)
        .def("CheckDofs", &VariableUtils::CheckDofs)
        .def("UpdateCurrentToInitialConfiguration", &VariableUtils::UpdateCurrentToInitialConfiguration)
        .def("UpdateInitialToCurrentConfiguration", &VariableUtils::UpdateInitialToCurrentConfiguration)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPosition)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariable)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariableAndPosition)
        ;


}

} // namespace Python.
} // Namespace Kratos
