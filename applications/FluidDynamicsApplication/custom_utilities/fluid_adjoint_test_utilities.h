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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAdjointTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using DofsVectorType =  std::vector<Dof<double>::Pointer>;

    using EquationIdVectorType = std::vector<std::size_t>;

    ///@}
    ///@name Static Operations
    ///@{

    template<class TDataType>
    static TDataType CalculateRelaxedVariableRate(
        const double BossakAlpha,
        const Variable<TDataType>& rVariable,
        const NodeType& rNode);

    ///@}
    ///@name Classes
    ///@{

    template<class TContainerType>
    class Testing
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

        template<class TDataType>
        static void RunAdjointEntityDerivativesTest(
            ModelPart& rPrimalModelPart,
            ModelPart& rAdjointModelPart,
            const std::function<void(ModelPart&)>& rUpdateModelPart,
            const Variable<TDataType>& rVariable,
            const std::function<void(Matrix&, typename TContainerType::data_type&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
            const IndexType EquationOffset,
            const IndexType DerivativeOffset,
            const double Delta,
            const double Tolerance);

        ///@}
    private:
        ///@name Private Operations
        ///@{

        static TContainerType& GetContainer(ModelPart& rModelPart);

        static void CalculateResidual(
            Vector& residual,
            typename TContainerType::data_type& rEntity,
            const ProcessInfo& rProcessInfo);

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
private:
    ///@name Private Operations
    ///@{

    template<class TDataType>
    static std::function<double&(NodeType&, const IndexType)> GetPerturbationMethod(
        const Variable<TDataType>& rPerturbationVariable);

    template<class TDataType>
    static IndexType GetVariableDimension(
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_TEST_UTILITIES_H_INCLUDED