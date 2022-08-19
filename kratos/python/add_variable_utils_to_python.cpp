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


/**
 * @brief Auxiliary set variable export function
 * This function is required to export the SetVariable overloaded method with a unique name
 * @tparam TDataType The variable data type
 * @tparam Variable<TDataType> The variable type
 * @param rVariableUtils Reference to the self variable utils class
 * @param rVariable Reference to the variable to be set
 * @param rValue Reference to the value to set
 * @param rNodes Reference to the nodes container
 */
template<class TDataType, class TVarType = Variable<TDataType>>
void VariableUtilsSetVariable(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    ModelPart::NodesContainerType &rNodes)
{
    rVariableUtils.SetVariable(rVariable, rValue, rNodes);
}

/**
 * @brief Auxiliary set variable export function
 * This function is required to export the SetVariable with flag overloaded method with a unique name
 * @tparam TDataType The variable data type
 * @tparam Variable<TDataType> The variable type
 * @param rVariableUtils Reference to the self variable utils class
 * @param rVariable Reference to the variable to be set
 * @param rValue Reference to the value to set
 * @param rNodes Reference to the nodes container
 * @param Flag Flag to filter the nodes that are set
 * @param CheckValue Flag value to be checked
 */
template <class TDataType, class TVarType = Variable<TDataType>>
void VariableUtilsSetVariableForFlag(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    ModelPart::NodesContainerType &rNodes,
    const Flags Flag,
    const bool CheckValue = true)
{
    rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);
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
#define KRATOS_VARIABLE_UTILS_TYPES bool, int, double, array_1d<double, 3>, array_1d<double, 4>, array_1d<double, 6>, array_1d<double, 9>, Quaternion<double>, Vector, Matrix

// Unified interface
template<class TVariableType>
void AddVariableUtilsCommonInterfaceIterator(pybind11::module& m) {
    m.def("CopyModelPartNodalVar",                      VariableUtilsCopyModelPartNodalVar<Variable<TVariableType>>);
    m.def("CopyModelPartNodalVar",                      VariableUtilsCopyModelPartNodalVarWithDestination<Variable<TVariableType>>);
    m.def("CopyModelPartNodalVarToNonHistoricalVar",    CopyModelPartNodalVarToNonHistoricalVar<Variable<TVariableType>>);
    m.def("CopyModelPartNodalVarToNonHistoricalVar",    CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<TVariableType>>);
    m.def("CopyModelPartElementalVar",                  &VariableUtils::CopyModelPartElementalVar<Variable<TVariableType>>);
    m.def("SetVariable",                                VariableUtilsSetVariable<TVariableType>);
    m.def("SetVariable",                                [](VariableUtils& rVariableUtils, const Variable<TVariableType>& rVariable, const TVariableType& rValue, ModelPart::NodesContainerType& rNodes, const unsigned int Step){rVariableUtils.SetVariable(rVariable, rValue, rNodes, Step);});
    m.def("SetVariable",                                VariableUtilsSetVariableForFlag<TVariableType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariable<TVariableType, ModelPart::NodesContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariable<TVariableType, ModelPart::ElementsContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariable<TVariableType, ModelPart::ConditionsContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariableForFlag<TVariableType, ModelPart::NodesContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariableForFlag<TVariableType, ModelPart::ElementsContainerType>);
    m.def("SetNonHistoricalVariable",                   VariableUtilsSetNonHistoricalVariableForFlag<TVariableType, ModelPart::ConditionsContainerType>);
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
        .def("SetVectorVar", VariableUtilsSetVariable<array_1d<double,3>>)
        .def("SetVectorVar", [](VariableUtils& rVariableUtils, const Variable<array_1d<double, 3>>& rVariable, const array_1d<double, 3>& rValue, ModelPart::NodesContainerType& rNodes, const unsigned int Step){rVariableUtils.SetVariable(rVariable, rValue, rNodes, Step);})
        .def("SetVectorVar", VariableUtilsSetVariableForFlag<array_1d<double, 3>>)
        .def("SetScalarVar", VariableUtilsSetVariable<double>)
        .def("SetScalarVar", [](VariableUtils& rVariableUtils, const Variable<double>& rVariable, const double& rValue, ModelPart::NodesContainerType& rNodes, const unsigned int Step){rVariableUtils.SetVariable(rVariable, rValue, rNodes, Step);})
        .def("SetScalarVar", VariableUtilsSetVariableForFlag<double>)

        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVectorVar", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)

        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<int>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<double>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<array_1d<double, 3>>)

        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::NodesContainerType>)

        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ElementsContainerType>)

        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ConditionsContainerType>)

        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::MasterSlaveConstraintContainerType>)

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
        
        .def("SaveScalarVar", &VariableUtils::SaveVariable<double>)              // To be removed
        .def("SaveVectorVar", &VariableUtils::SaveVariable<array_1d<double, 3>>) // To be removed
        
        .def("SaveScalarNonHistoricalVar", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType>)              // To be removed
        .def("SaveVectorNonHistoricalVar", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>) // To be removed

        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        
        .def("CopyScalarVar", &VariableUtils::CopyVariable<double>)                                                                    // To be removed
        .def("CopyVectorVar", &VariableUtils::CopyVariable<array_1d<double, 3>>)                                                       // To be removed

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
