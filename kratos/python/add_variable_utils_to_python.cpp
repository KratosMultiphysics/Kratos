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
#include "python/add_variable_utils_to_python.h"
#include "includes/define_python.h"
#include "processes/process.h"

// Variable utilities
#include "utilities/variable_utils.h"

namespace Kratos {
namespace Python {

template<class TDataType>
void AddCopyModelPartFlaggedInterface(pybind11::module& m)
{
    m.def("CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool, const unsigned int,const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, ModelPart&, const Flags&, const bool, const unsigned int,const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedElementVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedElementVar));
    m.def("CopyModelPartFlaggedElementVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedElementVar));
    m.def("CopyModelPartFlaggedConditionVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedConditionVar));
    m.def("CopyModelPartFlaggedConditionVar", (void(VariableUtils::*)(const Variable<TDataType>&, const Variable<TDataType>&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedConditionVar));

    m.def("CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool, const unsigned int,const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, ModelPart&, const Flags&, const bool, const unsigned int))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar));
    m.def("CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar", (void(VariableUtils::*)(const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar));
    m.def("CopyModelPartFlaggedElementVar", (void(VariableUtils::*)(const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedElementVar));
    m.def("CopyModelPartFlaggedConditionVar", (void(VariableUtils::*)(const Variable<TDataType>&, const ModelPart&, ModelPart&, const Flags&, const bool))(&VariableUtils::CopyModelPartFlaggedConditionVar));
}

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

template< class TVarType >
void ApplyFixity(
    VariableUtils &rVariableUtils,
    const TVarType& rVar,
    const bool IsFixed,
    ModelPart::NodesContainerType& rNodes)
{
    rVariableUtils.ApplyFixity(rVar, IsFixed, rNodes);
}

template< class TVarType >
void ApplyFlaggedFixity(
    VariableUtils &rVariableUtils,
    const TVarType& rVar,
    const bool IsFixed,
    ModelPart::NodesContainerType& rNodes,
    const Flags& rFlag,
    const bool CheckValue)
{
    rVariableUtils.ApplyFixity(rVar, IsFixed, rNodes, rFlag, CheckValue);
}

template <class TDataType, class TContainerType, class TVarType = Variable<TDataType>>
void VariableUtilsSetNonHistoricalVariable(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    TContainerType &rContainer)
{
    rVariableUtils.SetNonHistoricalVariable(rVariable, rValue, rContainer);
}

template <class TDataType, class TContainerType, class TVarType = Variable<TDataType>>
void VariableUtilsSetNonHistoricalVariableForFlag(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    TContainerType &rContainer,
    const Flags Flag,
    const bool CheckValue = true)
{
    rVariableUtils.SetNonHistoricalVariable(rVariable, rValue, rContainer, Flag, CheckValue);
}

// List of valid interface variables
#define KRATOS_VARIABLE_UTILS_TYPES             bool, int, double, array_1d<double, 3>, array_1d<double, 4>, array_1d<double, 6>, array_1d<double, 9>, Quaternion<double>, Vector, Matrix
#define KRATOS_VARIABLE_UTILS_CONTAINERS        ModelPart::NodesContainerType, ModelPart::ElementsContainerType, ModelPart::ConditionsContainerType
#define KRATOS_VARIABLE_UTILS_CONTAINERS_EXT    ModelPart::MasterSlaveConstraintContainerType

// Unified interface
template<class TVariableType>
void AddVariableUtilsCommonInterfaceIterator(pybind11::module& m) {
    m.def("CopyModelPartNodalVar",                      [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const ModelPart &rOriginModelPart, ModelPart &rDestinationModelPart, const unsigned int BuffStep = 0){rVariableUtils.CopyModelPartNodalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);});
    m.def("CopyModelPartNodalVar",                      [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const Variable<TVariableType>& rDestinationVariable, const ModelPart &rOriginModelPart, ModelPart &rDestinationModelPart, const unsigned int BuffStep = 0){rVariableUtils.CopyModelPartNodalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);});
    m.def("CopyModelPartNodalVarToNonHistoricalVar",    [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const ModelPart &rOriginModelPart, ModelPart &rDestinationModelPart, const unsigned int BuffStep = 0){rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);});
    m.def("CopyModelPartNodalVarToNonHistoricalVar",    [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const Variable<TVariableType>& rDestinationVariable, const ModelPart &rOriginModelPart, ModelPart &rDestinationModelPart, const unsigned int BuffStep = 0){rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);});
    m.def("CopyModelPartElementalVar",                  &VariableUtils::CopyModelPartElementalVar<Variable<TVariableType>>);
    m.def("SetVariable",                                [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const TVariableType& rValue, ModelPart::NodesContainerType& rNodes){rVariableUtils.SetVariable(rVariable, rValue, rNodes);});
    m.def("SetVariable",                                [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const TVariableType& rValue, ModelPart::NodesContainerType& rNodes, const unsigned int Step){rVariableUtils.SetVariable(rVariable, rValue, rNodes, Step);});
    m.def("SetVariable",                                [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const TVariableType &rValue, ModelPart::NodesContainerType& rNodes, const Flags Flag, const bool CheckValue = true){rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);});
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariable<TVariableType, ModelPart::NodesContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariable<TVariableType, ModelPart::ElementsContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariable<TVariableType, ModelPart::ConditionsContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariableForFlag<TVariableType, ModelPart::NodesContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariableForFlag<TVariableType, ModelPart::ElementsContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariableForFlag<TVariableType, ModelPart::ConditionsContainerType>);
    m.def("SetNonHistoricalVariableToZero",             &VariableUtils::SetNonHistoricalVariableToZero<TVariableType, ModelPart::NodesContainerType>);
    m.def("SetNonHistoricalVariableToZero",             &VariableUtils::SetNonHistoricalVariableToZero<TVariableType, ModelPart::ElementsContainerType>);
    m.def("SetNonHistoricalVariableToZero",             &VariableUtils::SetNonHistoricalVariableToZero<TVariableType, ModelPart::ConditionsContainerType>);
    m.def("SetNonHistoricalVariableToZero",             &VariableUtils::SetNonHistoricalVariableToZero<TVariableType, ModelPart::MasterSlaveConstraintContainerType>);
    m.def("SetHistoricalVariableToZero",                &VariableUtils::SetHistoricalVariableToZero<TVariableType>);
    m.def("SaveVariable",                               &VariableUtils::SaveVariable<TVariableType>);
    m.def("SaveNonHistoricalVariable",                  &VariableUtils::SaveNonHistoricalVariable<TVariableType, ModelPart::NodesContainerType>);
    m.def("SaveNonHistoricalVariable",                  &VariableUtils::SaveNonHistoricalVariable<TVariableType, ModelPart::ElementsContainerType>);
    m.def("SaveNonHistoricalVariable",                  &VariableUtils::SaveNonHistoricalVariable<TVariableType, ModelPart::ConditionsContainerType>);
    m.def("CopyVariable",                               &VariableUtils::CopyVariable<TVariableType>);

    AddCopyModelPartFlaggedInterface<TVariableType>(m);
}

template<typename... args>
void AddVariableUtilsCommonInterface(pybind11::module& m) {
    (AddVariableUtilsCommonInterfaceIterator<args>(m), ...);
}

void AddVariableUtilsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    auto python_variable_utils = py::class_<VariableUtils>(m, "VariableUtils")
        .def(py::init<>())
        
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::NodesContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ConditionsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ElementsContainerType>)

        .def("WeightedAccumulateConditionVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ConditionsContainerType, int>)
        .def("WeightedAccumulateConditionVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ConditionsContainerType, int>)
        .def("WeightedAccumulateElementVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ElementsContainerType, int>)
        .def("WeightedAccumulateElementVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ElementsContainerType, int>)
        .def("WeightedAccumulateConditionVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ConditionsContainerType, double>)
        .def("WeightedAccumulateConditionVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ConditionsContainerType, double>)
        .def("WeightedAccumulateElementVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ElementsContainerType, double>)
        .def("WeightedAccumulateElementVariableOnNodes", &VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ElementsContainerType, double>)
        
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
        
        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        
        .def("ApplyFixity", ApplyFixity<Variable<double>>)
        .def("ApplyFixity", ApplyFlaggedFixity<Variable<double>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<Variable<double>>)

        .def("GetSolutionStepValuesVector", py::overload_cast<
                            const ModelPart::NodesContainerType&,
                            const Variable<array_1d<double,3>>&,
                            const unsigned int,
                            const unsigned int>(&VariableUtils::GetSolutionStepValuesVector))
        .def("GetSolutionStepValuesVector", py::overload_cast<
                            const ModelPart::NodesContainerType&,
                            const Variable<double>&,
                            const unsigned int>(&VariableUtils::GetSolutionStepValuesVector))
        .def("SetSolutionStepValuesVector", py::overload_cast<
                            ModelPart::NodesContainerType&,
                            const Variable<array_1d<double,3>>&,
                            const Vector&,
                            const unsigned int>(&VariableUtils::SetSolutionStepValuesVector))
        .def("SetSolutionStepValuesVector", py::overload_cast<
                            ModelPart::NodesContainerType&,
                            const Variable<double>&,
                            const Vector&,
                            const unsigned int>(&VariableUtils::SetSolutionStepValuesVector))

        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalVariable<double>)
        .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalVariable<array_1d<double, 3>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumNonHistoricalNodeVectorVariable", &VariableUtils::SumNonHistoricalNodeVectorVariable)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<Variable<double>>)
        .def("SumConditionVectorVariable", &VariableUtils::SumConditionVectorVariable)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<Variable<double>>)
        .def("SumElementVectorVariable", &VariableUtils::SumElementVectorVariable)
        .def("AddDof", &VariableUtils::AddDof<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<Variable<double>>)
        .def_static("AddDofsList", &VariableUtils::AddDofsList)
        .def_static("AddDofsList", &VariableUtils::AddDofsWithReactionsList)
        .def("CheckVariableKeys", &VariableUtils::CheckVariableKeys)
        .def("UpdateCurrentToInitialConfiguration", &VariableUtils::UpdateCurrentToInitialConfiguration)
        .def("UpdateInitialToCurrentConfiguration", &VariableUtils::UpdateInitialToCurrentConfiguration)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPosition)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariable)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariableAndPosition)
        .def("GetCurrentPositionsVector", &VariableUtils::GetCurrentPositionsVector)
        .def("GetInitialPositionsVector", &VariableUtils::GetInitialPositionsVector)
        .def("SetCurrentPositionsVector", &VariableUtils::SetCurrentPositionsVector)
        .def("SetInitialPositionsVector", &VariableUtils::SetInitialPositionsVector)
        ;

    AddVariableUtilsCommonInterface<KRATOS_VARIABLE_UTILS_TYPES>(m);
}

} // namespace Python.
} // Namespace Kratos
