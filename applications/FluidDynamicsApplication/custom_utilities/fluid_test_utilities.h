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

#pragma once

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes

namespace Kratos
{
///@name Classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using PropertiesType = ModelPart::PropertiesType;

    using DofsVectorType =  std::vector<Dof<double>::Pointer>;

    using EquationIdVectorType = std::vector<std::size_t>;

    using NodesContainerType = ModelPart::NodesContainerType;

    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    using ElementsContainerType = ModelPart::ElementsContainerType;

    ///@}
    ///@name Static operations
    ///@{

    template<class TDataType>
    static void AssignRandomValues(
        TDataType& rValue,
        const std::string& rSeed,
        const int DomainSize,
        const double MinValue = 0.0,
        const double MaxValue = 1.0);

    static ModelPart& CreateTestModelPart(
        Model& rModel,
        const std::string& rModelPartName,
        const std::string& rElementName,
        const std::string& rConditionName,
        const std::function<void(PropertiesType&)>& rSetElementProperties,
        const std::function<void(PropertiesType&)>& rSetConditionProperties,
        const std::function<void(ModelPart&)>& rAddNodalSolutionStepVariablesFuncion,
        const std::function<void(NodeType&)>& rAddDofsFunction,
        const int BufferSize = 2);

    template<class TDataType>
    static void RandomFillHistoricalVariable(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const double MinValue = 0.0,
        const double MaxValue = 1.0,
        const int Step = 0)
    {
        RandomFillHistoricalVariable<TDataType>(rModelPart, rVariable, rVariable.Name(), MinValue, MaxValue, Step);
    }

    template<class TDataType>
    static void RandomFillHistoricalVariable(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const std::string& rSeedExtension,
        const double MinValue = 0.0,
        const double MaxValue = 1.0,
        const int Step = 0);

    template<class TContainerType, class TDataType>
    static void RandomFillNonHistoricalVariable(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize,
        const double MinValue = 0.0,
        const double MaxValue = 1.0)
    {
        RandomFillNonHistoricalVariable<TContainerType, TDataType>(rContainer, rVariable, rVariable.Name(), DomainSize, MinValue, MaxValue);
    }

    template<class TContainerType, class TDataType>
    static void RandomFillNonHistoricalVariable(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const std::string& rSeedExtension,
        const IndexType DomainSize,
        const double MinValue = 0.0,
        const double MaxValue = 1.0);

    template<class TContainerType>
    static void RunEntityGetDofListTest(
        const TContainerType& rContainer,
        const ProcessInfo& rProcessInfo,
        const std::vector<const Variable<double>*>& rDofVariablesList);

    template<class TContainerType>
    static void RunEntityEquationIdVectorTest(
        const TContainerType& rContainer,
        const ProcessInfo& rProcessInfo,
        const std::vector<const Variable<double>*>& rDofVariablesList);

    template<class TContainerType>
    static void RunEntityGetValuesVectorTest(
        const TContainerType& rContainer,
        const std::vector<const Variable<double>*>& rDofVariablesList);

    template<class TContainerType>
    static void RunEntityGetFirstDerivativesVectorTest(
        const TContainerType& rContainer,
        const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod);

    template<class TContainerType>
    static void RunEntityGetSecondDerivativesVectorTest(
        const TContainerType& rContainer,
        const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod);

    ///@}
    ///@name Classes
    ///@{

    template<class TContainerType>
    class GetContainer
    {
    public:
        TContainerType& operator()(ModelPart& rModelPart);
    };

    ///@}
};

///@}

} // namespace Kratos