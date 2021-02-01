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

#if !defined(KRATOS_FLUID_ADJOINT_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_FLUID_ADJOINT_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>
#include <vector>

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

class FluidAdjointTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using DofsVectorType =  std::vector<Dof<double>::Pointer>;

    using EquationIdVectorType = std::vector<std::size_t>;

    ///@}
    ///@name Classes
    ///@{

    template <class TDataType>
    class DataTypeUtilities
    {
    public:
        ///@name Static Operations
        ///@{

        static TDataType ComputeRelaxedVariableRate(
            const double BossakAlpha,
            const Variable<TDataType>& rVariable,
            const NodeType& rNode);

        static std::function<double&(NodeType&, const IndexType)> GetPerturbationMethod(
            const Variable<TDataType>& rPerturbationVariable);

        static IndexType GetVariableDimension(
            const Variable<TDataType>& rVariable,
            const ProcessInfo& rProcessInfo);

        static void AssignRandomValues(
            TDataType& rValue,
            const std::string& rSeed,
            const IndexType DomainSize,
            const double MinValue = 0.0,
            const double MaxValue = 1.0);

        ///@}
    };

    template<class TContainerType, class TDataType>
    class ContainerDataTypeUtilities
    {
    public:
        ///@name Static Operations
        ///@{

        static void RunAdjointEntityTest(
            ModelPart& rPrimalModelPart,
            ModelPart& rAdjointModelPart,
            const std::function<void(ModelPart&)>& rUpdateModelPart,
            const Variable<TDataType>& rVariable,
            const std::function<void(Matrix&, typename TContainerType::data_type&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
            const IndexType EquationOffset,
            const IndexType DerivativeOffset,
            const double Delta,
            const double Tolerance);

        static void RandomFillNodalHistoricalVariable(
            ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const double MinValue = 0.0,
            const double MaxValue = 1.0,
            const int Step = 0);

        static void RandomFillContainerVariable(
            ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const double MinValue = 0.0,
            const double MaxValue = 1.0);

        ///@}
    };

    template<class TContainerType>
    class ContainerUtilities
    {
    public:
        ///@name Static Operations
        ///@{

        static void RunAdjointEntityGetDofListTest(
            ModelPart& rModelPart,
            const std::vector<const Variable<double>*>& rDofVariablesList);

        static void RunAdjointEntityEquationIdVectorTest(
            ModelPart& rModelPart,
            const std::vector<const Variable<double>*>& rDofVariablesList);

        static void RunAdjointEntityGetValuesVectorTest(
            ModelPart& rModelPart,
            const std::vector<const Variable<double>*>& rDofVariablesList);

        static void RunAdjointEntityGetFirstDerivativesVectorTest(
            ModelPart& rModelPart,
            const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod);

        static void RunAdjointEntityGetSecondDerivativesVectorTest(
            ModelPart& rModelPart,
            const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod);

        static TContainerType& GetContainer(ModelPart& rModelPart);

        ///@}
    };

    ///@}
    ///@name Static operations
    ///@{

    template <class TClassType>
    static void CalculateResidual(
        Vector& residual,
        TClassType& rClassTypeObject,
        const ProcessInfo& rProcessInfo);

    static ModelPart& CreateTestModelPart(
        Model& rModel,
        const std::string& rElementName,
        const std::string& rConditionName,
        const std::function<ModelPart::PropertiesType::Pointer(ModelPart&)>& rGetElementProperties,
        const std::function<ModelPart::PropertiesType::Pointer(ModelPart&)>& rGetConditionProperties,
        const std::function<void(ModelPart&)>& rAddNodalSolutionStepVariablesFuncion,
        const std::function<void(NodeType&)>& rAddDofsFunction,
        const int BufferSize = 2);

    template<class TDataType>
    static void RandomFillNodalHistoricalVariable(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const double MinValue = 0.0,
        const double MaxValue = 1.0,
        const int Step = 0)
    {
        ContainerDataTypeUtilities<ModelPart::ElementsContainerType, TDataType>::RandomFillNodalHistoricalVariable(
            rModelPart, rVariable, MinValue, MaxValue, Step);
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_TEST_UTILITIES_H_INCLUDED